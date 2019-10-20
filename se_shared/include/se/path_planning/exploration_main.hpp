//
// Created by anna on 02/10/19.
//

#ifndef SUPEREIGHT_EXPLORATION_MAIN_HPP
#define SUPEREIGHT_EXPLORATION_MAIN_HPP

#include "candidate_view.hpp"
#include "path_planner_exploration.hpp"
#include "se/planner_config.h"

namespace se {

namespace exploration {
/**
 * calculates exploration path
 * for developing and visualization purposes also return candidate views
 * @tparam T Ofusion or SDF
 * @param shared_ptr octree_ptr : increments reference count and makes the function an owner, the
 * octree will stay alive as long as the function is sing it
 * @param frontier_map
 * @param res [m/voxel]
 * @param step interval to stride along a ray
 * @param planning_config
 * @param config supereight configuration
 * @param pose current camera pose
 * @param [out] path
 * @param [out] cand_views
 */
template<typename T>
int getExplorationPath(std::shared_ptr<Octree<T>> octree_ptr,
                       const float res,
                       const float step,
                       const Planning_Configuration &planning_config,
                       const Configuration &config,
                       const Eigen::Matrix4f &pose,
                       const Eigen::Vector3i &lower_bound,
                       const Eigen::Vector3i &upper_bound,
                       const float ground_height,
                       set3i &frontier_map,
                       VecPose &path,
                       VecPose &cand_views) {

  pose3D start = getCurrPose(pose, res);
  int valid_path = -1;
  auto collision_checker_v =
      aligned_shared < CollisionCheckerV < T > > (octree_ptr, planning_config);

  LOG(INFO) << "frontier map size " << frontier_map.size();
  // Candidate view generation
  CandidateView<T> candidate_view
      (octree_ptr, planning_config, collision_checker_v, res, config, pose, step, ground_height);
  std::shared_ptr<CandidateView<T>> cand_view_ptr= std::make_shared<CandidateView<T>>
      (*candidate_view);
  int has_candidate_positions =
      candidate_view.getCandidateViews(frontier_map, planning_config.frontier_cluster_size);
    // skip the path planning part
    PathPlannerExploration<T> path_planner_exploration(octree_ptr,
                                                       planning_config,
                                                       res,
                                                       config,
                                                       pose,
                                                       ground_height,
                                                       lower_bound,
                                                       upper_bound,
                                                       cand_view_ptr);
  valid_path = path_planner_exploration.planPathToCandidates();
  LOG(INFO)<< "path planned " << valid_path;

  candidate_view.calculateCandidateViewGain();
 LOG(INFO)<< "calc cand view "  ;
  int best_cand_idx = -1;
  bool use_curr_pose = true;
  bool force_travelling = candidate_view.curr_pose_.information_gain < candidate_view.getTargetIG();
  if (valid_path>0) {
    best_cand_idx = candidate_view.getBestCandidate();
    DLOG(INFO) << "[se/candview] best candidate is "
               << (candidate_view.candidates_[best_cand_idx].pose.p * res).format(InLine);
    // std::cout << " path length of best cand "
    //           << candidate_view.candidates_[best_cand_idx].path.size() << std::endl;
    use_curr_pose =
        candidate_view.candidates_[best_cand_idx].utility < candidate_view.curr_pose_.utility;
    LOG(INFO) << "force travelling " << force_travelling << " use curr pose " << use_curr_pose
              << " candidate utility " << candidate_view.candidates_[best_cand_idx].utility
              << " curr pose utility " << candidate_view.curr_pose_.utility;

  }

  VecPose path_tmp;
  if (valid_path>0 && (!use_curr_pose || force_travelling)) {
    // candidate_view.addPathSegments(planning_config.robot_safety_radius*2.5 ,best_cand_idx);
    // path_tmp = candidate_view.candidates_[best_cand_idx].path;
    path_tmp = path_planner_exploration.getFinalPath(candidate_view.candidates_[best_cand_idx]);
    // best_candidate = candidate_view.candidates_[best_cand_idx];
  } else {
    path_tmp = path_planner_exploration.getFinalPath(candidate_view.curr_pose_);
    // best_candidate = candidate_view.curr_pose_;
  }
  for (int i = 0; i <= planning_config.num_cand_views; i++) {
    if (candidate_view.candidates_[i].pose.p == Eigen::Vector3f(0, 0, 0)) {
      continue;
    }
    cand_views.push_back(candidate_view.candidates_[i].pose);
  }

  for (const auto &pose : path_tmp) {
    DLOG(INFO) << "cand view " << (pose.p * res).format(InLine) << " " << pose.q.w() << " "
               << pose.q.vec().format(InLine);
    path.push_back(pose);
  }
  if (candidate_view.getExplorationStatus() == 1 || frontier_map.size() == 0) {
    return 1;
  } else {
    return -1;
  }

}
} // exploration
} // se

#endif //SUPEREIGHT_EXPLORATION_MAIN_HPP
