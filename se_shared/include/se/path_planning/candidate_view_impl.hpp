

#ifndef SUPEREIGHT_CANDIDATE_VIEW_IMPL_HPP
#define SUPEREIGHT_CANDIDATE_VIEW_IMPL_HPP

#include "candidate_view.hpp"
namespace se {

namespace exploration {

/**
 * helper function to print the face voxels
 * TODO sparsify
 * @param _voxel center voxel cooridnates
 */
template<typename T>
void CandidateView<T>::printFaceVoxels(const Eigen::Vector3i &_voxel) {
  VecVec3i face_neighbour_voxel(6);
  face_neighbour_voxel[0] << _voxel.x() - 1, _voxel.y(), _voxel.z();
  face_neighbour_voxel[1] << _voxel.x() + 1, _voxel.y(), _voxel.z();
  face_neighbour_voxel[2] << _voxel.x(), _voxel.y() - 1, _voxel.z();
  face_neighbour_voxel[3] << _voxel.x(), _voxel.y() + 1, _voxel.z();
  face_neighbour_voxel[4] << _voxel.x(), _voxel.y(), _voxel.z() - 1;
  face_neighbour_voxel[5] << _voxel.x(), _voxel.y(), _voxel.z() + 1;
  std::cout << "[se/cand view] face voxel states ";
  for (const auto &face_voxel : face_neighbour_voxel) {
    std::cout << octree_ptr_->get(face_voxel).st << " ";
  }
  std::cout << std::endl;

}
template<typename T>
int CandidateView<T>::getCandidateViews(set3i &frontier_blocks_map,
                                        const int frontier_cluster_size) {
  int counter = 0;
  while (getNumValidCandidates() == 0 && frontier_blocks_map.size() != 0) {
    DLOG(INFO) << "get candidates";
    sampleCandidateViews(frontier_blocks_map, frontier_cluster_size);

    if (counter == planning_config_.num_cand_views) {
      LOG(INFO) << "no candidates ";
      return -1;
    }
    counter++;
  }
  return 1;
}
/**
 *
 * generates random valid candidate views from the frontier map
 * @tparam T
 * @param frontier_blocks_map
 */
template<typename T>
void CandidateView<T>::sampleCandidateViews(set3i &frontier_blocks_map,
                                            const int frontier_cluster_size) {

  mapvec3i frontier_voxels_map;
  node_iterator<T> node_it(*octree_ptr_);

  // get all frontier voxels inside a voxel block
  for (auto it = frontier_blocks_map.begin(); it != frontier_blocks_map.end(); it++) {
    VecVec3i frontier_voxels = node_it.getFrontierVoxels(*it);

    if (!frontier_voxels.empty()) {
      frontier_voxels_map[*it] = frontier_voxels;
    }
  }

  if (frontier_voxels_map.empty()) {
    candidates_.clear();
    LOG(INFO) << "No frontier voxels left. Exploration done.";
    return;
  }

  int sampling_step = std::ceil(frontier_voxels_map.size() / num_sampling_);
  if (sampling_step == 0) sampling_step = 1;
  DLOG(INFO) << "sampling every " << sampling_step << "th morton code";

  std::random_device rd;
  std::default_random_engine generator(planning_config_.random_generator_seed);
  if (planning_config_.random_generator_seed == 0) {
    std::default_random_engine generator(rd());
  }

  std::uniform_int_distribution<int> distribution_block(0, frontier_voxels_map.size() - 1);
  auto it = frontier_voxels_map.begin();
#pragma omp parallel for
  for (int i = 0; i < num_sampling_; i++) {

    std::advance(it, sampling_step);
    uint64_t rand_morton = it->first;

    if (frontier_voxels_map[rand_morton].size() < frontier_cluster_size
        || frontier_voxels_map[rand_morton].empty()) {
      continue;
    }

    std::uniform_int_distribution<int>
        distribution_voxel(0, frontier_voxels_map[rand_morton].size() - 1);
    // random frontier voxel inside the the randomly chosen voxel block
    int rand_voxel = distribution_voxel(generator);
    DLOG(INFO) << " rand voxel " << rand_voxel << " size "
               << frontier_voxels_map[rand_morton].size();
    Eigen::Vector3i candidate_frontier_voxel = frontier_voxels_map[rand_morton].at(rand_voxel);
    boundHeight(candidate_frontier_voxel.z(),
                planning_config_.planning_height_max + ground_height_,
                planning_config_.planning_height_min + ground_height_,
                res_);

    candidates_[i].pose.p = candidate_frontier_voxel.cast<float>();
    num_cands_++;
    DLOG(INFO) << "Cand voxel " << candidate_frontier_voxel.format(InLine);
  }
}

template<typename T>
float CandidateView<T>::getViewInformationGain(pose3D &pose) {
  LOG(INFO) << pose.p.format(InLine);
  float gain = 0.0;
  const float r_max = planning_config_.sensor_depth; //from camera model
  const float fov_hor = planning_config_.fov_hor;
  const float fov_vert = planning_config_.fov_vert; // image size

  // temporary
  const int n_col = fov_vert / dphi_;
  const int n_row = 360 / dtheta_;

  float phi_rad, theta_rad;

  Eigen::MatrixXf gain_matrix(n_row, n_col + 2);
  Eigen::MatrixXf depth_matrix(n_row, n_col + 1);
  std::map<int, float> gain_per_yaw;
  gain_per_yaw.clear();
  Eigen::Vector3f vec(0.0, 0.0, 0.0); //[m]
  Eigen::Vector3f dir(0.0, 0.0, 0.0);// [m]
  Eigen::Vector3f cand_view_m = pose.p * res_;

  // sparse ray cast every x0 deg
  int row = 0;
  for (int theta = -180; theta < 180; theta += dtheta_) {
    theta_rad = static_cast<float>(M_PI * theta / 180.0); // deg to rad
    gain = 0.0;
    int col = 0;
    gain_matrix(row, col) = theta;
    depth_matrix(row, col) = theta;
    for (int phi = static_cast<int>(90 - fov_vert / 2); phi < 90 + fov_vert / 2; phi += dphi_) {
      col++;
      // calculate ray cast direction
      phi_rad = static_cast<float>(M_PI * phi / 180.0f);
      vec[0] = cand_view_m[0] + r_max * cos(theta_rad) * sin(phi_rad);
      vec[1] = cand_view_m[1] + r_max * sin(theta_rad) * sin(phi_rad);
      vec[2] = cand_view_m[2] + r_max * cos(phi_rad);
      dir = (vec - cand_view_m).normalized();
      // initialize ray
      se::ray_iterator<T> ray(*octree_ptr_, cand_view_m, dir, nearPlane, r_max);
      ray.next();
      // lower bound dist from camera
      const float t_min = ray.tcmin(); /* Get distance to the first intersected block */
      //get IG along ray
      std::pair<float, float> gain_tmp =
          t_min > 0.f ? getRayInformationGain(dir, t_min, r_max, cand_view_m) : std::make_pair(0.f,
                                                                                               0.f);
      gain_matrix(row, col) = gain_tmp.first;
      depth_matrix(row, col) = gain_tmp.second;
      gain += gain_tmp.first;
    }

    gain_matrix(row, col + 1) = gain;
    row++;
    // add gain to map
    if (!gain_per_yaw.insert(std::make_pair(theta, gain)).second) {
      LOG(INFO) << "Insertion failed. Key was present";
    }

  }

  int best_yaw = 0;
  float best_yaw_gain = 0.0;
  // best yaw evaluation
  getBestYawandGain(fov_hor, gain_per_yaw, best_yaw, best_yaw_gain);

  float yaw_rad = M_PI * best_yaw / 180.f;
  pose.q = toQuaternion(yaw_rad, 0.0, 0.0);
  return best_yaw_gain;
//    std::cout << "[se/candview] gain_matrix \n" << gain_matrix.transpose() << std::endl;
//    std::cout << "[se/candview] depth_matrix \n" << depth_matrix.transpose() << std::endl;
//    std::cout << "[se/candview] for cand " << cand_view.format(InLine) << " best theta angle is "
//              << best_yaw << ", best ig is " << best_yaw_gain << std::endl;

//    saveMatrixToDepthImage(depth_matrix.block(0, 1, n_row, n_col).transpose().rowwise().reverse(),
//                           cand_num,
//                           true);
//    saveMatrixToDepthImage(gain_matrix.block(0, 1, n_row, n_col).transpose().rowwise().reverse(),
//                           cand_num,
//                           false);
}

template<typename T>
void CandidateView<T>::getBestYawandGain(const float fov_hor,
                                         const std::map<int, float> &gain_per_yaw,
                                         int &best_yaw,
                                         float &best_yaw_gain) {

  for (int yaw = -180; yaw < 180; yaw++) {
    float yaw_score = 0.0;
// gain FoV horizontal
    for (int fov = -fov_hor / 2; fov < fov_hor / 2; fov++) {
      int theta = yaw + fov;
// wrap angle
      wrapYawDeg(theta);
      if (gain_per_yaw.count(theta)) {
        yaw_score += gain_per_yaw.at(theta);
      }
    }

    if (best_yaw_gain < yaw_score) {
      best_yaw_gain = yaw_score;
      best_yaw = yaw;
    }
  }
}

/**
 * march along the ray until max sensor depth or a surface is hit to calculate the entropy reduction
 *
 * @param direction [m]
 * @param tnear
 * @param tfar max sensor depth
 * @param origin
 * @return
 */
template<typename T>
std::pair<float, float> CandidateView<T>::getRayInformationGain(const Eigen::Vector3f &direction,
                                                                const float tnear,
                                                                const float tfar,
                                                                const Eigen::Vector3f &origin) {
  auto select_occupancy = [](typename VoxelBlock<T>::value_type val) { return val.x; };
  // march from camera away
  float ig_entropy = 0.0f;
  float t = tnear; // closer bound to camera
  if (tnear < tfar) {

    // occupancy prob in log2
    float prob_log = octree_ptr_->interp(origin + direction * t, select_occupancy);
    float weight = 1.0f;
    // check if current pos is free
    if (prob_log <= SURF_BOUNDARY + occ_thresh_) {
      ig_entropy = weight * getEntropy(prob_log);
      for (; t < tfar; t += step_) {
        const Eigen::Vector3f pos = origin + direction * t;
        const Eigen::Vector3i pos_v = (pos / res_).cast<int>();
        if (pos_v.x() <= 0 || pos_v.y() <= 0 || pos_v.z() <= 0 || pos_v.x() >= octree_ptr_->size()
            || pos_v.y() >= octree_ptr_->size() || pos_v.z() >= octree_ptr_->size())
          break;
        const auto data = octree_ptr_->get(pos_v);
        prob_log = octree_ptr_->interp(pos_v.cast<float>(), select_occupancy);
        ig_entropy += getEntropy(prob_log);

// next step along the ray hits a surface with a secure threshold return
        if (prob_log > SURF_BOUNDARY + occ_thresh_) {
          break;
        }
      }
    }
  } else {
    // TODO 0.4 is set fix in ray_iterator
    float num_it = (tfar - nearPlane) / step_;
    ig_entropy = num_it * 1.f;
    t = tfar;
  }

  return std::make_pair(ig_entropy, t);
}

// information gain calculation
// source [1] aeplanner gainCubature
template<typename T>
void CandidateView<T>::calculateCandidateViewGain() {
  float ig_sum = 0.f;
  int cand_counter = 0;
#pragma omp parallel for
  for (int i = 0; i < num_sampling_; i++) {

    if (candidates_[i].pose.p == Eigen::Vector3f(0, 0, 0)) {
      // cand_counter++;
      continue;
    }

    cand_counter++;
    float information_gain = getViewInformationGain(candidates_[i].pose);
    candidates_[i].information_gain = information_gain;
    LOG(INFO) << "cand " << i << " IG " << information_gain;
    ig_sum += candidates_[i].information_gain;
  }
  float information_gain = getViewInformationGain(curr_pose_.pose);
  curr_pose_.information_gain = information_gain;
  ig_sum += curr_pose_.information_gain;
  if (ig_sum / cand_counter < ig_target_) {
    exploration_status_ = 1;
    LOG(INFO) << "low information gain  " << ig_sum / cand_counter;

  }

}

template<typename T>
void CandidateView<T>::calculateUtility(Candidate &candidate) {
  const float max_yaw_rate = planning_config_.max_yaw_rate;
  float yaw_diff = toEulerAngles(candidate.pose.q).yaw - toEulerAngles(pose_.q).yaw;
  float t_path = candidate.path_length * res_;
  //wrap yaw
  wrapYawRad(yaw_diff);
  float t_yaw = fabs(yaw_diff / max_yaw_rate);
  candidate.utility = candidate.information_gain / (t_yaw + t_path);

  if (t_path == 0 && t_yaw < 0.001) {
    candidate.utility = 0;
  }

  DLOG(INFO) << "IG: " << candidate.information_gain << ", t_yaw: " << t_yaw << ", t_path: "
             << t_path << ", utility: " << candidate.utility;
}
/**
 * finds the best candidate based on information gain
 * @param cand_list
 * @return
 */
template<typename T>
int CandidateView<T>::getBestCandidate() {

  float max_utility = 0.0f;
  int best_cand_idx = 0;

  // path cost = voxel *res[m/vox] / (v =1m/s) = time
#pragma omp parallel for
  for (int i = 0; i < num_sampling_; i++) {
    if (candidates_[i].pose.p != Eigen::Vector3f(0, 0, 0)
        && candidates_[i].planning_solution_status >= 0) {
      calculateUtility(candidates_[i]);
    }
  }
  // highest utility in the beginning

  std::sort(candidates_.begin(),
            candidates_.end(),
            [](const auto &a, const auto &b) { return (a.utility > b.utility); });

  calculateUtility(curr_pose_);
  return 0;
}
} // exploration
} // se
#endif