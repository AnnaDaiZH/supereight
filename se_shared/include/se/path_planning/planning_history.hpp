/**
 * Information-theoretic exploration
 *
 * Copyright (C) 2019 Imperial College London.
 * Copyright (C) 2019 ETH Zürich.
 *
 * @file planning_history.hpp
 * @author Anna Dai
 * @date August 22, 2019
 */

#ifndef SUPEREIGHT_PLANNING_HISTORY_HPP
#define SUPEREIGHT_PLANNING_HISTORY_HPP

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


#include "se/octree.hpp"
#include "se/node_iterator.hpp"
#include "se/constant_parameters.h"
#include "se/utils/math_utils.h"
#include "se/config.h"
#include "se/planner_config.h"
#include "se/utils/eigen_utils.h"
#include "exploration_utils.hpp"


namespace se {

namespace exploration {

template<typename FieldType>
class PlanningHistoryManager {
 public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
   PlanningHistoryManager(){};
   PlanningHistoryManager(const std::shared_ptr<Octree<FieldType> > octree_ptr,
   	                      const Planning_Configuration &planning_config);


   void updateHistoryPath(const Candidate &path, const VecPose &path_short);
   void insertNewCandidates(const VecCandidate &candidates);
   void updateValidCandidates(const Eigen::Matrix4f &curr_pose);
   int useHistoryPath(const VecPose &path_short);

   VecCandidate& getOldCandidates(){return candidates_old_;}
   Candidate getLastPlannedTrajectory() {return last_path_;}
   VecCandidate getPathHistory(){return path_history_;}
   int getLocalMinimaCounter()const {return local_minima_counter_;}

   VecPose getOldPath();

 private:

  std::shared_ptr<Octree < FieldType> > octree_ptr_ = nullptr;
  Planning_Configuration planning_params_;
  Candidate last_path_;
  VecCandidate path_history_;
  VecCandidate candidates_old_;
  int local_minima_counter_;

  float voxel_dim_;
};


template<typename FieldType>
PlanningHistoryManager<FieldType>::PlanningHistoryManager(const std::shared_ptr<Octree <FieldType> > octree_ptr,
	const Planning_Configuration &planning_config)
:
octree_ptr_(octree_ptr),
planning_params_(planning_config),
local_minima_counter_(0){
	voxel_dim_ = octree_ptr->voxelDim();
}


template<typename FieldType>
void PlanningHistoryManager<FieldType>::updateHistoryPath(const Candidate &path, const VecPose &path_short){
  last_path_.clear();
  last_path_ = Candidate(path);
  path_history_.push_back(path);
  path_history_.at(path_history_.size()-1).path = path_short;
}

template<typename FieldType>
void PlanningHistoryManager<FieldType>::insertNewCandidates(const VecCandidate &candidates){
  LOG(INFO)<< candidates.size();
   if(candidates_old_.empty()){
   	for (const auto& new_cand : candidates)
   	{
   		candidates_old_.push_back(new_cand);
      DLOG(INFO) << candidates_old_.back().pose.p.format(InLine);
   	}
   	return;
   }
// check if candidate already exists if not. insert
   bool  is_new_cand = true;
	 for(const auto& new_cand : candidates ){
       for(auto old_cand : candidates_old_){
       	if(new_cand.pose.p == old_cand.pose.p || new_cand.pose.p == Eigen::Vector3f(0.f, 0.f,0.f)){
       		is_new_cand =false;
       	}
       }
       if(is_new_cand) candidates_old_.push_back(new_cand);
	}
}

template<typename FieldType>
void PlanningHistoryManager<FieldType>::updateValidCandidates(const Eigen::Matrix4f &curr_pose){
	Eigen::Vector3f curr_position = curr_pose.block<3,1>(0,3);
	Eigen::Vector3i curr_position_v = (curr_position/ voxel_dim_).cast<int>();
	if(candidates_old_.empty()) return;
	// check if the candidate is still a frontier
  DLOG(INFO) << candidates_old_.size();
  for( auto it = candidates_old_.begin() ; it != candidates_old_.end();){
    Eigen::Vector3f dist = curr_position_v.cast<float>() - it->pose.p;
    DLOG(INFO) << it->pose.p.format(InLine);
  	// if(octree_ptr_->get(it->pose.p.template cast<int>()).st != voxel_state::kFrontier
      // || dist.norm() < farPlane*2/voxel_dim_){
    if( dist.norm() < farPlane*2/voxel_dim_){
      DLOG(INFO) << it->pose.p.format(InLine);
      it = candidates_old_.erase(it);

  	} else{it++;}

  	}
DLOG(INFO) << candidates_old_.size();
}




/**
* if the robot is stuck in the same area for too long
* let it fly along the last path back
* A) to the position with the highest IG stored
* b) recalculate some IG => need CandidateView class
* determine max distance, num pose to evaluate
*/
template<typename FieldType>
int PlanningHistoryManager<FieldType>::useHistoryPath(const VecPose &path_short){
  Eigen::Vector3f dist_to_last_goal = path_short.back().p-path_history_.back().path.back().p;
// check if the end point is in the local minima area
  if(dist_to_last_goal.norm() < planning_params_.local_minima_radius/ voxel_dim_){
    if(local_minima_counter_>= 2){
  // if the counter exceeds threshold, send old path back
      return 1;
    } else{
  // if no, increase counter
      local_minima_counter_ ++;
      return -1;
    }
  } else{
    local_minima_counter_ =0;
    return -1;
  }

}


template<typename FieldType>
VecPose PlanningHistoryManager<FieldType>::getOldPath(){
  VecPose path_out;

  float max_utility=0.f;
  float curr_utility= 0.f;
  int cand_idx_from_back = 0;
  int iteration_considered = path_history_.size() < 6 ? path_history_.size() : 6;

  for(int i = 1 ; i <= iteration_considered ; i ++){
    curr_utility = path_history_.at(path_history_.size() -i).utility;
    // prevent to use another old path
    if(curr_utility == -1.f){
      break;
    }
    if(max_utility < curr_utility){
      max_utility = curr_utility;
      cand_idx_from_back = i;
    }
    std::cout << "max ut " << max_utility << " curr ut "<< curr_utility << " cand " << cand_idx_from_back << std::endl;
  }

  for(int i = 1; i<= cand_idx_from_back; i++){
    VecPose path_tmp = path_history_.at(path_history_.size() -i).path;
    for(const auto & pose : path_tmp){
      path_out.push_back(pose);
    }
  }

  Candidate tmp;
  tmp.path = path_out;
  tmp.information_gain = path_history_.at(path_history_.size() -cand_idx_from_back).information_gain;

  tmp.pose = path_history_.at(path_history_.size() -cand_idx_from_back).pose;
  path_history_.push_back(tmp);


  return path_out;
}

 } //exploration
} //se



#endif