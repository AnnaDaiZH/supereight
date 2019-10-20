/**
 * Information-theoretic exploration
 *
 * Copyright (C) 2019 Imperial College London.
 * Copyright (C) 2019 ETH Zürich.
 *
 * @file candidate_view.hpp
 * @author Anna Dai
 * @date August 22, 2019
 */

#ifndef SUPEREIGHT_CANDIDATE_VIEW_HPP
#define SUPEREIGHT_CANDIDATE_VIEW_HPP

#include <set>
#include <map>
#include <cstdlib>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <vector>
#include <random>
#include <iterator>
#include <type_traits>
#include <cmath>
#include <glog/logging.h>
#include <Eigen/StdVector>

#include "se/geometry/octree_collision.hpp"
#include "se/octree.hpp"
#include "se/node_iterator.hpp"
#include "se/constant_parameters.h"
#include "se/ray_iterator.hpp"
#include "se/utils/math_utils.h"
#include "se/config.h"
#include "se/utils/eigen_utils.h"
#include "exploration_utils.hpp"

#include "se/path_planner_ompl.hpp"


//typedef SE_FIELD_TYPE FieldType;
namespace se {

namespace exploration {

/**
 * Candidate View
 * all in voxel coord
 */



template<typename T>
class CandidateView {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef std::shared_ptr<CandidateView> SPtr;
  CandidateView(const std::shared_ptr<Octree<T> > octree_ptr,
                const Planning_Configuration &planning_config,
                const std::shared_ptr<CollisionCheckerV<T> > pcc,
                const float res,
                const Configuration &config,
                const Eigen::Matrix4f &curr_pose,
                const float step,
                const float ground_height);

  int getCandidateViews(set3i &frontier_blocks_map, const int frontier_cluster_size);
  void sampleCandidateViews(set3i &frontier_blocks_map, const int frontier_cluster_size);

  void printFaceVoxels(const Eigen::Vector3i &voxel);

  std::pair<float, float> getRayInformationGain(const Eigen::Vector3f &direction,
                                                const float tnear,
                                                const float tfar,
                                                const Eigen::Vector3f &origin);
  float getViewInformationGain(pose3D &pose);

  void calculateCandidateViewGain();

  void calculateUtility(Candidate &candidate);
  int getBestCandidate();

  void getBestYawandGain(const float fov_hor,
                         const std::map<int, float> &gain_per_yaw,
                         int &best_yaw,
                         float &best_yaw_gain);

  VecPose getFinalPath(Candidate &candidate);
  VecPose getYawPath(const pose3D &start, const pose3D &goal);
  VecPose fusePath(VecPose &path_tmp, VecPose &yaw_path);
  VecPose addPathSegments(const pose3D &start, const pose3D &goal);

  int getExplorationStatus() const { return exploration_status_; }

  float getTargetIG() const { return ig_target_; }
  float getTotalIG() const {return ig_total_; }
  int getNumValidCandidates() const { return num_cands_; }

  VecCandidate candidates_;
  Candidate curr_pose_;
 private:
  pose3D pose_;
  VecVec3i cand_views_;

  float res_; // [m/vx]
  Planning_Configuration planning_config_;
  Configuration config_;
  float occ_thresh_ = 0.f;
  float ig_total_;
  float ig_target_;
  int dtheta_;
  int dphi_;
  float step_;
  float ground_height_;
  int exploration_status_;

  int num_sampling_;
  int num_cands_;

  std::shared_ptr<CollisionCheckerV<T> > pcc_ = nullptr;
  std::shared_ptr<Octree<T> > octree_ptr_ = nullptr;
};

template<typename T>
CandidateView<T>::CandidateView(const std::shared_ptr<Octree<T> > octree_ptr,
                                const Planning_Configuration &planning_config,
                                const std::shared_ptr<CollisionCheckerV<T> > pcc,
                                const float res,
                                const Configuration &config,
                                const Eigen::Matrix4f &curr_pose,
                                const float step,
                                const float ground_height)
    :
    octree_ptr_(octree_ptr),
    planning_config_(planning_config),
    res_(res),
    config_(config),
    pcc_(pcc),
    dtheta_(planning_config.dtheta),
    dphi_(planning_config.dphi),
    step_(step),
    exploration_status_(0),
    num_sampling_(planning_config.num_cand_views),
    num_cands_(0),
    ground_height_(ground_height) {

  curr_pose_.pose = getCurrPose(curr_pose, res_);
  curr_pose_.path.push_back(curr_pose_.pose);
  curr_pose_.path_length = 0.f;
  pose_ = curr_pose_.pose;
  int n_col = planning_config.fov_vert / planning_config.dphi;
  int n_row = planning_config.fov_hor / planning_config.dtheta;
  ig_total_ = n_col * n_row * (farPlane / step) * getEntropy(0);
  ig_target_ = n_col * n_row * (farPlane / step) * getEntropy(log2(0.1 / (1.f - 0.1)));
  LOG(INFO)<< "ig total " << ig_total_ << " ig target " << ig_target_ ;
  candidates_.resize(num_sampling_);
  DLOG(INFO) << "setup candidate view";
}

} // namespace exploration
} // namespace se

#include "candidate_view_impl.hpp"
#endif //SUPEREIGHT_CANDIDATE_VIEW_HPP
