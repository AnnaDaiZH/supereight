#include <fstream>
#include <string>
#include <random>
#include "se/octree.hpp"
#include "se/node_iterator.hpp"
#include "se/node.hpp"
#include "se/utils/math_utils.h"
#include "se/utils/morton_utils.hpp"
#include "se/utils/eigen_utils.h"
#include "se/config.h"
#include "se/planner_config.h"
#include "se/path_planning/exploration_utils.hpp"
#include "se/path_planning/planning_history.hpp"
#include "se/utils/memory_pool.hpp"

#include "se/path_planner_ompl.hpp"
#include "se/ompl/collision_checker_voxel.hpp"
#include "se/ompl/motion_validator_occupancy_skeleton.hpp"

#include "gtest/gtest.h"

class CollisionUnitTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    tree_ = std::make_shared<se::Octree<OFusion> >();
    std::string
        filename = "/home/anna/exploration_ws/src/supereight/se_ompl/test/collision/w_box.bin";
    tree_->load(filename);
    std::cout << "file loaded" << std::endl;
    max_depth_ = log2(tree_->size());
    leaves_level_ = tree_->leaf_level();
    dim_ = tree_->voxelDim();
    planner_config_ = getDefaultPlanningConfig();

  }

  map3i createMap3i(se::MemoryPool<se::VoxelBlock<OFusion> > &block_buffer) {
    map3i morton_map;
    for (int i = 0; i < block_buffer.size(); i++) {
      const Eigen::Vector3i block_coord = block_buffer[i]->coordinates();
      const key_t morton_code = block_buffer[i]->code_;
      morton_map.emplace(morton_code, block_coord);
    }
    return morton_map;
  }

  std::shared_ptr<se::Octree<OFusion> > tree_;
  int max_depth_;
  int leaves_level_;
  int dim_;
  Planning_Configuration planner_config_;

};

TEST_F(CollisionUnitTest, CorridorCheckPass) {

  auto &block_buffer_base = tree_->getBlockBuffer();
  // std::shared_ptr<se::Octree<OFusion> > octree_ptr(&tree_);
  auto collision_checker =
      aligned_shared<se::exploration::CollisionCheckerV<OFusion> >(tree_, planner_config_);

  Eigen::Vector3i start = {89, 75, 72};
  Eigen::Vector3i end = {95, 75, 72};

  bool is_collision_free = collision_checker->isSegmentFlightCorridorSkeletonFree(start, end, 0, 3);
  EXPECT_TRUE(is_collision_free);

}

TEST_F(CollisionUnitTest, CorridorCheckFail) {

  auto &block_buffer_base = tree_->getBlockBuffer();
  // std::shared_ptr<se::Octree<OFusion> > octree_ptr(&tree_);
  auto collision_checker =
      aligned_shared<se::exploration::CollisionCheckerV<OFusion> >(tree_, planner_config_);

  Eigen::Vector3i start = {60, 50, 72};
  Eigen::Vector3i end = {80, 50, 72};

  bool is_collision_free = collision_checker->isSegmentFlightCorridorSkeletonFree(start, end, 0, 4);
  EXPECT_FALSE(is_collision_free);

}

TEST_F(CollisionUnitTest, PathPlanningPass) {

  auto &block_buffer_base = tree_->getBlockBuffer();
  // std::shared_ptr<se::Octree<OFusion> > octree_ptr(&tree_);
  auto collision_checker =
      aligned_shared<se::exploration::CollisionCheckerV<OFusion> >(tree_, planner_config_);
  auto path_planner_ompl_ptr = aligned_shared<se::exploration::PathPlannerOmpl<OFusion> >(tree_,
                                                                                          collision_checker,
                                                                                          planner_config_,
                                                                                          12.1f);
  map3i free_map = createMap3i(block_buffer_base);
  Eigen::Vector3i lower_bound ;
  Eigen::Vector3i upper_bound ;
  getFreeMapBounds(tree_, free_map, lower_bound, upper_bound);
  bool setup_planner = path_planner_ompl_ptr->setupPlanner(lower_bound, upper_bound);
  Eigen::Vector3i start = {89, 75, 72};
  Eigen::Vector3i end = {95, 75, 72};
  int solution_status = path_planner_ompl_ptr->planPath(start, end);
  EXPECT_EQ(solution_status, 1);

}

TEST_F(CollisionUnitTest, PathPlanningAroundWall) {

  auto &block_buffer_base = tree_->getBlockBuffer();
  // std::shared_ptr<se::Octree<OFusion> > octree_ptr(&tree_);
  auto collision_checker =
      aligned_shared<se::exploration::CollisionCheckerV<OFusion> >(tree_, planner_config_);
  auto path_planner_ompl_ptr = aligned_shared<se::exploration::PathPlannerOmpl<OFusion> >(tree_,
                                                                                          collision_checker,
                                                                                          planner_config_,
                                                                                          12.1f);
  map3i free_map = createMap3i(block_buffer_base);
  Eigen::Vector3i lower_bound ;
  Eigen::Vector3i upper_bound ;
  getFreeMapBounds(tree_, free_map, lower_bound, upper_bound);
  bool setup_planner = path_planner_ompl_ptr->setupPlanner(lower_bound, upper_bound);
  Eigen::Vector3i start = {80, 80, 72};
  Eigen::Vector3i end = {90, 50, 72};
  int solution_status = path_planner_ompl_ptr->planPath(start, end);
  EXPECT_EQ(solution_status, 1);

}


TEST_F(CollisionUnitTest, HistoryCandidateUpdate) {

  auto &block_buffer_base = tree_->getBlockBuffer();
  se::exploration::PlanningHistoryManager<OFusion> history_manager(tree_, planner_config_);
  VecCandidate old_candidates;
  VecCandidate new_candidates;
  Candidate tmp;
  tmp.pose.p = Eigen::Vector3f(80, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(83, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(85, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(81, 80, 72);
  new_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(0, 0, 0);
  new_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(83, 80, 72);
  new_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(0, 0, 0);
  new_candidates.push_back(tmp);


  history_manager.insertNewCandidates(old_candidates);
  EXPECT_EQ(history_manager.getOldCandidates().size(), 3 );

  history_manager.insertNewCandidates(new_candidates);

 EXPECT_EQ(history_manager.getOldCandidates().size(), 4);

}

TEST_F(CollisionUnitTest, UpdateHistoryPath){

  se::exploration::PlanningHistoryManager<OFusion> history_manager(tree_, planner_config_);
  VecPose path;
  path.push_back({ {0.f, 0.f, 0.f}, {1.f, 0.f, 0.f, 0.f}});
  path.push_back({ {0.2f, 0.2f, 0.2f}, {1.f, 0.f, 0.f, 0.f}});
  path.push_back({ {0.3f, 0.2f, 0.2f}, {1.f, 0.f, 0.f, 0.f}});

  Candidate best_candidate;
  best_candidate.path.push_back({ {0.f, 0.f, 0.f}, {1.f, 0.f, 0.f, 0.f}});
  best_candidate.path.push_back({ {0.2f, 0.2f, 0.2f}, {1.f, 0.f, 0.f, 0.f}});
  best_candidate.path.push_back({ {0.3f, 0.2f, 0.2f}, {1.f, 0.f, 0.f, 0.f}});
  best_candidate.path.push_back({ {0.4f, 0.2f, 0.2f}, {1.f, 0.f, 0.f, 0.f}});

  history_manager.updateHistoryPath(best_candidate, path);

  EXPECT_EQ(history_manager.getLastPlannedTrajectory().path.size(), best_candidate.path.size());
  EXPECT_EQ(history_manager.getPathHistory().back().path.size(), path.size());
}

TEST_F(CollisionUnitTest, UpdateValidCandidatesDeleteAll) {


  se::exploration::PlanningHistoryManager<OFusion> history_manager(tree_, planner_config_);

  Eigen::Matrix4f curr_pose;
  curr_pose << 1 , 0,0,1 ,0,1,0,1, 0,0,1,1, 0,0,0,1;
  VecCandidate old_candidates;

  Candidate tmp;
  tmp.pose.p = Eigen::Vector3f(80, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(83, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(85, 80, 72);
  old_candidates.push_back(tmp);

  history_manager.updateValidCandidates(curr_pose);

  EXPECT_EQ(history_manager.getOldCandidates().size(), 0 );
  history_manager.insertNewCandidates(old_candidates);
  EXPECT_EQ(history_manager.getOldCandidates().size(), 3 );
  history_manager.updateValidCandidates(curr_pose);
  EXPECT_EQ(history_manager.getOldCandidates().size(), 0);
}

TEST_F(CollisionUnitTest, UpdateValidCandidates2Frontiers) {
  voxel_traits<OFusion>::value_type data;
  data.st = voxel_state::kFrontier;
  tree_->set(80,80,72, data);
  tree_->set(85,80,72, data);

  se::exploration::PlanningHistoryManager<OFusion> history_manager(tree_, planner_config_);

  Eigen::Matrix4f curr_pose;
  curr_pose << 1 , 0,0,1 ,0,1,0,1, 0,0,1,1, 0,0,0,1;
  VecCandidate old_candidates;

  Candidate tmp;
  tmp.pose.p = Eigen::Vector3f(80, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(83, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(85, 80, 72);
  old_candidates.push_back(tmp);

  history_manager.updateValidCandidates(curr_pose);

  EXPECT_EQ(history_manager.getOldCandidates().size(), 0 );
  history_manager.insertNewCandidates(old_candidates);
  EXPECT_EQ(history_manager.getOldCandidates().size(), 3 );
  history_manager.updateValidCandidates(curr_pose);
  EXPECT_EQ(history_manager.getOldCandidates().size(), 2);
}

TEST_F(CollisionUnitTest, UpdateValidCandidatesNoMapUpdate) {
  voxel_traits<OFusion>::value_type data;
  data.st = voxel_state::kFrontier;
  tree_->set(80,80,72, data);
  tree_->set(85,80,72, data);

  se::exploration::PlanningHistoryManager<OFusion> history_manager(tree_, planner_config_);

  Eigen::Matrix4f curr_pose;
  curr_pose << 1 , 0,0,1 ,0,1,0,1, 0,0,1,1, 0,0,0,1;
  VecCandidate old_candidates;

  Candidate tmp;
  tmp.pose.p = Eigen::Vector3f(80, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(83, 80, 72);
  old_candidates.push_back(tmp);
  tmp.pose.p = Eigen::Vector3f(85, 80, 72);
  old_candidates.push_back(tmp);

  history_manager.updateValidCandidates(curr_pose);

  EXPECT_EQ(history_manager.getOldCandidates().size(), 0 );
  history_manager.insertNewCandidates(old_candidates);
  EXPECT_EQ(history_manager.getOldCandidates().size(), 3 );
  history_manager.updateValidCandidates(curr_pose);
  EXPECT_EQ(history_manager.getOldCandidates().size(), 2);
}