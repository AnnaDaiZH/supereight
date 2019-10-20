//
// Created by anna on 02/10/19.
//

#ifndef EXPLORATION_WS_PATH_PLANNER_EXPLORATION_IMPL_HPP
#define EXPLORATION_WS_PATH_PLANNER_EXPLORATION_IMPL_HPP

#include "path_planner_exploration.hpp"

namespace se {

namespace exploration {

template<typename T>
PathPlannerExploration<T>::PathPlannerExploration(const std::shared_ptr<Octree < T>>
octree_ptr,
const Planning_Configuration &planning_config,
const float res,
const Configuration &config,
const Eigen::Matrix4f &curr_pose,
const float ground_height,
const Eigen::Vector3i &lower_bound,
const Eigen::Vector3i &upper_bound,
std::shared_ptr<CandidateView < T>> candidate_view_ptr)
:

octree_ptr_ (octree_ptr), planning_config_(planning_config), res_(res), config_(config), ground_height_(
    ground_height), lower_bound_(lower_bound), upper_bound_(upper_bound), candidate_view_ptr_(
    candidate_view_ptr) {

  pose_ = getCurrPose(curr_pose, res_);
  DLOG(INFO) << "setup path planner exploration";

};

template<typename T>
int PathPlannerExploration<T>::planPathToCandidates() {

  int path_planned = -1;
#pragma omp parallel for
  for (int i = 0; i < planning_config_.num_cand_views; i++) {
    LOG(INFO) << "cand num " << i;
    if (candidate_view_ptr_->candidates_[i].pose.p != Eigen::Vector3f(0, 0, 0)) {
      auto collision_checker = aligned_shared<CollisionCheckerV<T> >(octree_ptr_, planning_config_);
      auto path_planner_ompl_ptr = aligned_shared<PathPlannerOmpl<T> >(octree_ptr_,
                                                                       collision_checker,
                                                                       planning_config_,
                                                                       ground_height_);
      // LOG(INFO) << "Candidate " << i << " goal coord " << cand_views[i].p.format(InLine);
      DLOG(INFO) << "Candidate " << i << " start " << pose_.p.format(InLine) << " goal "
                 << candidate_view_ptr_->candidates_[i].pose.p.format(InLine);

      bool setup_planner = path_planner_ompl_ptr->setupPlanner(lower_bound_, upper_bound_);
      DLOG(INFO) << "setup planner successful " << setup_planner;
      path_planned = path_planner_ompl_ptr->planPath(pose_.p.cast<int>(),
                                                     candidate_view_ptr_->candidates_[i].pose.p.template cast<
                                                         int>());
      candidate_view_ptr_->candidates_[i].planning_solution_status = path_planned;
      DLOG(INFO) << "path planned " << path_planned;

      if (path_planned < 0) {
        Planning_Configuration planning_config_tmp = planning_config_;
        float robot_safety_radius = planning_config_.robot_safety_radius;
        for (; robot_safety_radius >= 0.3f; robot_safety_radius -= 0.2f) {
          planning_config_tmp.robot_safety_radius = robot_safety_radius;
          DLOG(INFO) << "robot_safety_radius " << planning_config_tmp.robot_safety_radius;
          path_planner_ompl_ptr = aligned_shared<PathPlannerOmpl<T> >(octree_ptr_,
                                                                      collision_checker,
                                                                      planning_config_tmp,
                                                                      ground_height_);

          setup_planner = path_planner_ompl_ptr->setupPlanner(lower_bound_, upper_bound_);

          path_planned = path_planner_ompl_ptr->planPath(pose_.p.cast<int>(),
                                                         candidate_view_ptr_->candidates_[i].pose.p.template cast<
                                                             int>());
          candidate_view_ptr_->candidates_[i].planning_solution_status = path_planned;

          if (path_planned > 0) {
            LOG(INFO) << "Found a path for cand " << i;
            break;
          } else {
            candidate_view_ptr_->candidates_[i].path_length =
                (pose_.p - candidate_view_ptr_->candidates_[i].pose.p).squaredNorm();
            VecPose vec;
            vec.push_back(candidate_view_ptr_->candidates_[i].pose);
            candidate_view_ptr_->candidates_[i].path = vec;
          }
        }
        // voxel coord

      }

      if (path_planned > 0) {
        candidate_view_ptr_->candidates_[i].path_length = path_planner_ompl_ptr->getPathLength();

        VecPose vec = path_planner_ompl_ptr->getPathSegments_m(); //[voxel coord]
        candidate_view_ptr_->candidates_[i].path = vec;

        if (path_planned == 2) {
          DLOG(INFO) << "cand changed from "
                     << candidate_view_ptr_->candidates_[i].pose.p.format(InLine) << " to "
                     << candidate_view_ptr_->candidates_[i].path[
                         candidate_view_ptr_->candidates_[i].path.size() - 1].p.format(InLine);
          candidate_view_ptr_->candidates_[i].pose.p = candidate_view_ptr_->candidates_[i].path[
              candidate_view_ptr_->candidates_[i].path.size() - 1].p;
        }
        // DLOG(INFO) << "seg length " << all_path[i].size() << std::endl;
      }
    }
  }

  return path_planned;
}

template<typename T>
VecPose PathPlannerExploration<T>::getFinalPath(Candidate &candidate) {

  VecPose path, path_tmp, path_yaw;
// first add points between paths
  if (candidate.path.size() > 2) {
    for (auto i = 1; i < candidate.path.size(); i++) {

      DLOG(INFO) << "segment start " << candidate.path[i - 1].p.format(InLine) << " end "
                 << candidate.path[i].p.format(InLine);
      // interpolate long path segments
      if ((candidate.path[i].p - candidate.path[i - 1].p).norm()
          > planning_config_.v_max * planning_config_.dt) {
        path_tmp = addPathSegments(candidate.path[i - 1], candidate.path[i]);
      } else {
        path_tmp.push_back(candidate.path[i - 1]);
        path_tmp.push_back(candidate.path[i]);
      }

      if (planning_config_.yaw_optimization || i == candidate.path.size() - 1) {
        candidate_view_ptr_->getViewInformationGain(candidate.path[i]);
      } else {
        // look in the flight direction
        float x_side = (candidate.path[i].p.x() - candidate.path[i - 1].p.x());
        float y_side = (candidate.path[i].p.y() - candidate.path[i - 1].p.y());
        auto angle = (float) atan2(y_side, x_side);
        wrapYawRad(angle);
        candidate.path[i].q = toQuaternion(angle, 0, 0);
        path_yaw = getYawPath(candidate.path[i - 1], candidate.path[i]);
      }

      if (i == 1) {
        path_yaw = getYawPath(pose_, candidate.path[i]);
      } else {
        path_yaw = getYawPath(candidate.path[i - 1], candidate.path[i]);
      }
      VecPose path_fused = fusePath(path_tmp, path_yaw);
      // push back the new path
      for (const auto &pose : path_fused) {
        path.push_back(pose);
      }
    }
  } else {
    path_yaw = getYawPath(pose_, candidate.pose);
    if ((candidate.pose.p - pose_.p).norm() > planning_config_.v_max * planning_config_.dt) {
      path_tmp = addPathSegments(pose_, candidate.pose);
    } else {
      path_tmp = candidate.path;
    }
    path = fusePath(path_tmp, path_yaw);
  }

  return path;
}

template<typename T>
VecPose PathPlannerExploration<T>::fusePath(VecPose &path_tmp, VecPose &yaw_path) {
  VecPose path_out;
  if (yaw_path.size() >= path_tmp.size()) {

    for (int i = 0; i < yaw_path.size(); i++) {
      if (i < path_tmp.size()) {
        yaw_path[i].p = path_tmp[i].p;
      } else {
        yaw_path[i].p = path_tmp[path_tmp.size() - 1].p;
      }
    }

    path_out = yaw_path;
  } else {

    for (int i = 0; i < path_tmp.size(); i++) {
      if (i < yaw_path.size()) {
        path_tmp[i].q = yaw_path[i].q;
      } else {
        DLOG(INFO) << "last yaw "
                   << toEulerAngles(yaw_path[yaw_path.size() - 1].q).yaw * 180 / M_PI;
        path_tmp[i].q = yaw_path[yaw_path.size() - 1].q;
      }
    }
    path_out = path_tmp;
  }
  return path_out;
}

template<typename T>
VecPose PathPlannerExploration<T>::getYawPath(const pose3D &start, const pose3D &goal) {
  VecPose path;
  pose3D pose_tmp;
  const float max_yaw_rate = planning_config_.max_yaw_rate;
//  path.push_back(start);
  float yaw_diff = toEulerAngles(goal.q).yaw - toEulerAngles(start.q).yaw;
  // LOG(INFO) << "[interpolateYaw] yaw diff " << yaw_diff ;
  wrapYawRad(yaw_diff);
  // LOG(INFO) << "[interpolateYaw] yaw diff " << yaw_diff ;
  // interpolate yaw
  float yaw_increment = (planning_config_.max_yaw_rate * planning_config_.dt) / std::abs(yaw_diff);
  for (float i = 0.f; i <= 1.0f; i += yaw_increment) {
    pose_tmp = goal;
    pose_tmp.q = start.q.slerp(i, goal.q);
    DLOG(INFO) << "intermediate yaw " << toEulerAngles(pose_tmp.q).yaw * 180 / M_PI;
    pose_tmp.q.normalize();
    path.push_back(pose_tmp);
  }

  path.push_back(goal);
  return path;
}

template<typename T>
VecPose PathPlannerExploration<T>::addPathSegments(const pose3D &start_in, const pose3D &goal_in) {

  VecPose path_out;
  pose3D start = start_in;
  pose3D goal = goal_in;
  int start_h = start.p.z();
  int goal_h = goal.p.z();
  boundHeight(start_h,
              planning_config_.planning_height_max + ground_height_,
              planning_config_.planning_height_min + ground_height_,
              res_);
  boundHeight(goal_h,
              planning_config_.planning_height_max + ground_height_,
              planning_config_.planning_height_min + ground_height_,
              res_);

  DLOG(INFO) << start.p.format(InLine) << " " << start_h;
  start.p.z() = (float) start_h;
  goal.p.z() = (float) goal_h;
  path_out.push_back(start);
  float dist = (goal.p - start.p).norm();
  float dist_increment = (planning_config_.v_max * planning_config_.dt / res_) / dist; // [voxel]
  Eigen::Vector3f dir = (goal.p - start.p).normalized();
  for (float t = 0.f; t <= 1.0f; t += dist_increment) {
    Eigen::Vector3f intermediate_point = start.p + dir * t * dist;
    pose3D tmp(intermediate_point, {1.f, 0.f, 0.f, 0.f});
    DLOG(INFO) << "intermediate_point " << intermediate_point.format(InLine);
    path_out.push_back(tmp);

  }
  path_out.push_back(goal);
  return path_out;

}

}
}
#endif //EXPLORATION_WS_PATH_PLANNER_EXPLORATION_IMPL_HPP
