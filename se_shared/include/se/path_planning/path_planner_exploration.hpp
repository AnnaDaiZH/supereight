//
// Created by anna on 02/10/19.
//

#ifndef EXPLORATION_WS_PATH_PLANNER_EXPLORATION_HPP
#define EXPLORATION_WS_PATH_PLANNER_EXPLORATION_HPP

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
#include "candidate_view.hpp"

namespace se {

namespace exploration {

template<typename T>
class PathPlannerExploration {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW typedef std::shared_ptr<PathPlannerExploration>
  SPtr;
  PathPlannerExploration(const std::shared_ptr<Octree<T>> octree_ptr,
                         const Planning_Configuration &planning_config,
                         const float res,
                         const Configuration &config,
                         const Eigen::Matrix4f &curr_pose,
                         const float ground_height,
                         const Eigen::Vector3i &lower_bound,
                         const Eigen::Vector3i &upper_bound,
                         std::shared_ptr<CandidateView<T>> candidate_view_ptr);

  int planPathToCandidates();
  VecPose getFinalPath(Candidate &candidate);
  VecPose getYawPath(const pose3D &start, const pose3D &goal);
  VecPose fusePath(VecPose &path_tmp, VecPose &yaw_path);
  VecPose addPathSegments(const pose3D &start, const pose3D &goal);

 private:
  pose3D pose_;
  std::shared_ptr<CandidateView<T>> candidate_view_ptr_= nullptr;
  float res_; // [m/vx]
  Planning_Configuration planning_config_;
  Configuration config_;
  float ground_height_;
  Eigen::Vector3i lower_bound_;
  Eigen::Vector3i upper_bound_;

  std::shared_ptr<Octree<T>> octree_ptr_ = nullptr;
};

} //exploration

}//se

#include "path_planner_exploration_impl.hpp"
#endif //EXPLORATION_WS_PATH_PLANNER_EXPLORATION_HPP
