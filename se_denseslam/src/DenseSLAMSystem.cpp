/*

 Copyright (c) 2014 University of Edinburgh, Imperial College, University of Manchester.
 Developed in the PAMELA project, EPSRC Programme Grant EP/K008730/1

 This code is licensed under the MIT License.


 Copyright 2016 Emanuele Vespa, Imperial College London

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <se/DenseSLAMSystem.h>
#include <se/ray_iterator.hpp>
#include <se/algorithms/meshing.hpp>
#include <se/geometry/octree_collision.hpp>
#include <se/vtk-io.h>
#include "preprocessing.cpp"
#include "tracking.cpp"
#include "rendering.cpp"
#include "bfusion/mapping_impl.hpp"
#include "kfusion/mapping_impl.hpp"
#include "bfusion/alloc_impl.hpp"
#include "kfusion/alloc_impl.hpp"
//#include "se/boundary_extraction.hpp"

#include "se/ompl/prob_collision_checker.hpp"
#include "se/utils/planning_parameter.hpp"
#include "se/path_planner_ompl.hpp"
PerfStats Stats;
static bool print_kernel_timing = false;

DenseSLAMSystem::DenseSLAMSystem(const Eigen::Vector2i &inputSize,
                                 const Eigen::Vector3i &volumeResolution,
                                 const Eigen::Vector3f &volumeDimensions,
                                 const Eigen::Vector3f &initPose,
                                 std::vector<int> &pyramid,
                                 const Configuration &config,
                                 const Planning_Configuration &planning_config)
    :
    DenseSLAMSystem(inputSize,
                    volumeResolution,
                    volumeDimensions,
                    se::math::toMatrix4f(initPose),
                    pyramid,
                    config,
                    planning_config) {}

DenseSLAMSystem::DenseSLAMSystem(const Eigen::Vector2i &inputSize,
                                 const Eigen::Vector3i &volumeResolution,
                                 const Eigen::Vector3f &volumeDimensions,
                                 const Eigen::Matrix4f &initPose,
                                 std::vector<int> &pyramid,
                                 const Configuration &config,
                                 const Planning_Configuration &planning_config)
    :
    computation_size_(inputSize),
    vertex_(computation_size_.x(), computation_size_.y()),
    normal_(computation_size_.x(), computation_size_.y()),
    float_depth_(computation_size_.x(), computation_size_.y()) {
  planning_config_ = planning_config;
  this->init_pose_ = initPose.block<3, 1>(0, 3);
  this->volume_dimension_ = volumeDimensions;
  this->volume_resolution_ = volumeResolution;
  this->mu_ = config.mu;
  pose_ = initPose;
  Tbc_ << 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1;
  raycast_pose_ = initPose * Tbc_;

  this->iterations_.clear();
  for (std::vector<int>::iterator it = pyramid.begin(); it != pyramid.end(); it++) {
    this->iterations_.push_back(*it);
  }

  viewPose_ = &pose_;

  if (getenv("KERNEL_TIMINGS"))
    print_kernel_timing = true;

  // internal buffers to initialize
  reduction_output_.resize(sizeof(reduction_output_.data()) * sizeof(reduction_output_.data()) * 4);
  tracking_result_.resize(computation_size_.x() * computation_size_.y());

  for (unsigned int i = 0; i < iterations_.size(); ++i) {
    int downsample = 1 << i;
    scaled_depth_.push_back(se::Image<float>(computation_size_.x() / downsample,
                                             computation_size_.y() / downsample));

    input_vertex_.push_back(se::Image<Eigen::Vector3f>(computation_size_.x() / downsample,
                                                       computation_size_.y() / downsample));

    input_normal_.push_back(se::Image<Eigen::Vector3f>(computation_size_.x() / downsample,
                                                       computation_size_.y() / downsample));
  }

  // ********* BEGIN : Generate the gaussian *************
  size_t gaussianS = radius * 2 + 1;
  gaussian_.reserve(gaussianS);
  int x;
  for (unsigned int i = 0; i < gaussianS; i++) {
    x = i - 2;
    gaussian_[i] = expf(-(x * x) / (2 * delta * delta));
  }

  // ********* END : Generate the gaussian *************

  discrete_vol_ptr_ = aligned_shared<se::Octree<FieldType> >();
  discrete_vol_ptr_->init(volume_resolution_.x(), volume_dimension_.x());
  volume_ =
      Volume<FieldType>(volume_resolution_.x(), volume_dimension_.x(), discrete_vol_ptr_.get());

}

bool DenseSLAMSystem::preprocessing(const unsigned short *inputDepth,
                                    const Eigen::Vector2i &inputSize,
                                    const bool filterInput) {

  mm2metersKernel(float_depth_, inputDepth, inputSize);
  if (filterInput) {
    bilateralFilterKernel(scaled_depth_[0], float_depth_, gaussian_, e_delta, radius);
  } else {
    std::memcpy(scaled_depth_[0].data(),
                float_depth_.data(),
                sizeof(float) * computation_size_.x() * computation_size_.y());
  }
  return true;
}

bool DenseSLAMSystem::tracking(const Eigen::Vector4f &k,
                               float icp_threshold,
                               unsigned tracking_rate,
                               unsigned frame) {

  if (frame % tracking_rate != 0)
    return false;

  // half sample the input depth maps into the pyramid levels
  for (unsigned int i = 1; i < iterations_.size(); ++i) {
    halfSampleRobustImageKernel(scaled_depth_[i], scaled_depth_[i - 1], e_delta * 3, 1);
  }

  // prepare the 3D information from the input depth maps
  Eigen::Vector2i localimagesize = computation_size_;
  for (unsigned int i = 0; i < iterations_.size(); ++i) {
    Eigen::Matrix4f invK = getInverseCameraMatrix(k / float(1 << i));
    depth2vertexKernel(input_vertex_[i], scaled_depth_[i], invK);
    if (k.y() < 0)
      vertex2normalKernel<true>(input_normal_[i], input_vertex_[i]);
    else
      vertex2normalKernel<false>(input_normal_[i], input_vertex_[i]);
    localimagesize /= 2;;
  }

  old_pose_ = pose_;
  const Eigen::Matrix4f projectReference = getCameraMatrix(k) * raycast_pose_.inverse();

  for (int level = iterations_.size() - 1; level >= 0; --level) {
    Eigen::Vector2i localimagesize
        (computation_size_.x() / (int) pow(2, level), computation_size_.y() / (int) pow(2, level));
    for (int i = 0; i < iterations_[level]; ++i) {

      trackKernel(tracking_result_.data(),
                  input_vertex_[level],
                  input_normal_[level],
                  vertex_,
                  normal_,
                  pose_,
                  projectReference,
                  dist_threshold,
                  normal_threshold);

      reduceKernel(reduction_output_.data(),
                   tracking_result_.data(),
                   computation_size_,
                   localimagesize);

      if (updatePoseKernel(pose_, reduction_output_.data(), icp_threshold))
        break;

    }
  }
  return checkPoseKernel(pose_,
                         old_pose_,
                         reduction_output_.data(),
                         computation_size_,
                         track_threshold);
}

bool DenseSLAMSystem::raycasting(const Eigen::Vector4f &k, float mu, unsigned int frame) {

  bool doRaycast = false;

  if (frame > 2) {
    raycast_pose_ = pose_ * Tbc_; // camera to world
    float step = volume_dimension_.x() / volume_resolution_.x();

    //TODO maintainability
    raycastKernel(volume_,
                  vertex_,
                  normal_,
                  raycast_pose_ * getInverseCameraMatrix(k),
                  nearPlane,
                  farPlane,
                  mu,
                  step,
                  step * BLOCK_SIDE);
//                  surface_voxel_set_,
//                  frontier_voxel_set_,
//                  occlusion_voxel_set_);
    doRaycast = true;
  }
  return doRaycast;
}

bool DenseSLAMSystem::integration(const Eigen::Vector4f &k,
                                  unsigned int integration_rate,
                                  float mu,
                                  unsigned int frame) {

  if (((frame % integration_rate) == 0) || (frame <= 3)) {

    float voxelsize = volume_._extent / volume_._size;
    int num_vox_per_pix = volume_._extent / ((se::VoxelBlock<FieldType>::side) * voxelsize);
    size_t total = num_vox_per_pix * computation_size_.x() * computation_size_.y();
    allocation_list_.reserve(total);

    unsigned int allocated = 0;
    if (std::is_same<FieldType, SDF>::value) {
      allocated = buildAllocationList(allocation_list_.data(),
                                      allocation_list_.capacity(),
                                      *volume_._map_index,
                                      pose_ * Tbc_,
                                      getCameraMatrix(k),
                                      float_depth_.data(),
                                      computation_size_,
                                      volume_._size,
                                      voxelsize,
                                      2 * mu);
    } else if (std::is_same<FieldType, OFusion>::value) {
      allocated = buildOctantList(allocation_list_.data(),
                                  allocation_list_.capacity(),
                                  *volume_._map_index,
                                  pose_,
                                  getCameraMatrix(k),
                                  float_depth_.data(),
                                  computation_size_,
                                  voxelsize,
                                  compute_stepsize,
                                  step_to_depth,
                                  mu);
    }

    // Allocate the nodes determined through the raycasting process
    volume_._map_index->allocate(allocation_list_.data(), allocated);

    if (std::is_same<FieldType, SDF>::value) {
      struct sdf_update funct(float_depth_.data(),
                              Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
                              mu,
                              100);
      se::functor::projective_map(*volume_._map_index,
                                  Sophus::SE3f(pose_).inverse(),
                                  getCameraMatrix(k),
                                  Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
                                  funct);
    } else if (std::is_same<FieldType, OFusion>::value) {

      float timestamp = (1.f / 30.f) * frame;

      struct bfusion_update funct(float_depth_.data(),
                                  Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
                                  mu,
                                  timestamp,
                                  voxelsize);
// Update all active nodes and voxels using the bfusion_update function

      se::functor::projective_map(*volume_._map_index,
                                  Sophus::SE3f(pose_).inverse(),
                                  getCameraMatrix(k),
                                  Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
                                  funct);
    }

    // if(frame % 15 == 0) {
    //   std::stringstream f;
    //   f << "./slices/integration_" << frame << ".vtk";
    //   save3DSlice(*volume_._map_index, Eigen::Vector3i(0, 200, 0),
    //       Eigen::Vector3i(volume_._size, 201, volume_._size),
    //       Eigen::Vector3i::Constant(volume_._size), f.str().c_str());
    //   f.str("");
    //   f.clear();
    // }
  } else {
    return false;
  }
  return true;
}
//
//bool DenseSLAMSystem::integration(const Eigen::Vector4f& k, unsigned int integration_rate,
//    float mu, unsigned int frame,
//    std::vector<Eigen::Vector3i> *occupied_voxels,
//    std::vector<Eigen::Vector3i> *freed_voxels) {
//
//  if (((frame % integration_rate) == 0) || (frame <= 3)) {
//
//    float voxelsize =  volume_._extent/volume_._size;
//    int num_vox_per_pix = volume_._extent/((se::VoxelBlock<FieldType>::side)*voxelsize);
//    size_t total = num_vox_per_pix * computation_size_.x() *
//      computation_size_.y();
//    allocation_list_.reserve(total);
//
//    float float_depth[computation_size_.x() * computation_size_.y()];
//    for (int y = 0; y < computation_size_.y(); y++) {
//      for (int x = 0; x < computation_size_.x(); x++) {
//        float_depth[x + y*computation_size_.x()] = float_depth_.data()[x + y*computation_size_.x()];
//      }
//    }
//
//    unsigned int allocated = 0;
//    if(std::is_same<FieldType, SDF>::value) {
//     allocated  = buildAllocationList(allocation_list_.data(),
//         allocation_list_.capacity(),
//        *volume_._map_index, pose_, getCameraMatrix(k), float_depth,
//        computation_size_, volume_._size,
//      voxelsize, 2*mu);
//    } else if(std::is_same<FieldType, OFusion>::value) {
//     allocated = buildOctantList(allocation_list_.data(), allocation_list_.capacity(),
//         *volume_._map_index,
//         pose_, getCameraMatrix(k), float_depth, computation_size_, voxelsize,
//         compute_stepsize, step_to_depth, 6*mu);
//    }
//
//    volume_._map_index->allocate(allocation_list_.data(), allocated);
//
//    if(std::is_same<FieldType, SDF>::value) {
//      struct sdf_update funct(float_depth,
//          Eigen::Vector2i(computation_size_.x(), computation_size_.y()), mu, 100);
//      se::functor::projective_map(*volume_._map_index,
//          Sophus::SE3f(pose_).inverse(),
//          getCameraMatrix(k),
//          Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
//          funct);
//    } else if(std::is_same<FieldType, OFusion>::value) {
//
//      float timestamp = (1.f/30.f)*frame;
//
//      // Initialize bfusion_update function to update each node/voxel updating step
//      struct bfusion_update funct(float_depth,
//          Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
//          mu, timestamp, voxelsize, occupied_voxels, freed_voxels);
//
//      // Update all active nodes and voxels using the bfusion_update function
//      se::functor::projective_map(*volume_._map_index,
//          Sophus::SE3f(pose_).inverse(),
//          getCameraMatrix(k),
//          Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
//          funct);
//    }
//
//    // if(frame % 15 == 0) {
//    //   std::stringstream f;
//    //   f << "./slices/integration_" << frame << ".vtk";
//    //   save3DSlice(*volume_._map_index, Eigen::Vector3i(0, 200, 0),
//    //       Eigen::Vector3i(volume_._size, 201, volume_._size),
//    //       Eigen::Vector3i::Constant(volume_._size), f.str().c_str());
//    //   f.str("");
//    //   f.clear();
//    // }
//  } else {
//    return false;
//  }
//  return true;
//}

// TODO fix alignment
bool DenseSLAMSystem::integration(const Eigen::Vector4f &k,
                                  unsigned int integration_rate,
                                  float mu,
                                  unsigned int frame,
                                  set3i *updated_blocks,
                                  set3i *frontier_blocks,
                                  set3i *free_blocks) {

  if (((frame % integration_rate) == 0) || (frame <= 3)) {

    float voxelsize = volume_._extent / volume_._size;
    int num_vox_per_pix = volume_._extent / ((se::VoxelBlock<FieldType>::side) * voxelsize);

    size_t total = num_vox_per_pix * computation_size_.x() * computation_size_.y();
    allocation_list_.reserve(total);

    float float_depth[computation_size_.x() * computation_size_.y()];
    for (int y = 0; y < computation_size_.y(); y++) {
      for (int x = 0; x < computation_size_.x(); x++) {
        float_depth[x + y * computation_size_.x()] =
            float_depth_.data()[x + y * computation_size_.x()];
      }
    }

    unsigned int allocated = 0;
    if (std::is_same<FieldType, SDF>::value) {
      allocated = buildAllocationList(allocation_list_.data(),
                                      allocation_list_.capacity(),
                                      *volume_._map_index,
                                      pose_ * Tbc_,
                                      getCameraMatrix(k),
                                      float_depth,
                                      computation_size_,
                                      volume_._size,
                                      voxelsize,
                                      2 * mu);
    } else if (std::is_same<FieldType, OFusion>::value) {
      allocated = buildOctantList(allocation_list_.data(),
                                  allocation_list_.capacity(),
                                  *volume_._map_index,
                                  pose_ * Tbc_,
                                  getCameraMatrix(k),
                                  float_depth,
                                  computation_size_,
                                  voxelsize,
                                  compute_stepsize,
                                  step_to_depth,
                                  6 * mu);
    }

    volume_._map_index->allocate(allocation_list_.data(), allocated);

    if (std::is_same<FieldType, SDF>::value) {
      struct sdf_update funct
          (float_depth, Eigen::Vector2i(computation_size_.x(), computation_size_.y()), mu, 100);
      se::functor::projective_map(*volume_._map_index,
                                  Sophus::SE3f(pose_ * Tbc_).inverse(),
                                  getCameraMatrix(k),
                                  Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
                                  funct);
    } else if (std::is_same<FieldType, OFusion>::value) {

      float timestamp = (1.f / 30.f) * frame;
      struct bfusion_update funct(float_depth_.data(),
                                  Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
                                  mu,
                                  timestamp,
                                  voxelsize,
                                  updated_blocks,
                                  free_blocks,
                                  frontier_blocks,
                                  false,
                                  1);
// Update all active nodes and voxels using the bfusion_update function

      se::functor::projective_map(*volume_._map_index,
                                  Sophus::SE3f(pose_ * Tbc_).inverse(),
                                  getCameraMatrix(k),
                                  Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
                                  funct);
      // struct bfusion_update frontier_funct(float_depth,
      //                                      Eigen::Vector2i(computation_size_.x(),
      //                                                      computation_size_.y()),
      //                                      mu,
      //                                      timestamp,
      //                                      voxelsize,
      //                                      frontier_blocks,
      //                                      true);

      // se::functor::projective_map(*volume_._map_index,
      //                             Sophus::SE3f(pose_ * Tbc_).inverse(),
      //                             getCameraMatrix(k),
      //                             Eigen::Vector2i(computation_size_.x(), computation_size_.y()),
      //                             frontier_funct);
      std::set<uint64_t> *copy_frontier_blocks = frontier_blocks;
      bool update_frontier_map = (frame % integration_rate) == 0;
      updateFrontierMap(volume_, frontier_map_, copy_frontier_blocks, update_frontier_map);
      insertBlocksToMap(free_map_, free_blocks);
      // std::cout << "[se/denseslam] free_map_  size  " << free_map_.size() << std::endl;
      std::cout << "[se/denseslam] frontier_map_ size " << frontier_map_.size() << std::endl;
    }

    // if(frame % 15 == 0) {
    //   std::stringstream f;
    //   f << "./slices/integration_" << frame << ".vtk";
    //   save3DSlice(*volume_._map_index, Eigen::Vector3i(0, 200, 0),
    //       Eigen::Vector3i(volume_._size, 201, volume_._size),
    //       Eigen::Vector3i::Constant(volume_._size), f.str().c_str());
    //   f.str("");
    //   f.clear();
    // }
  } else {
    return false;
  }
  return true;
}

bool DenseSLAMSystem::planning(VecPose &path,
                               VecPose &cand_views,
                               mapvec3i *free_blocks,
                               int *exploration_done) {
  se::exploration::initNewPosition(pose_ * Tbc_,
                                   planning_config_,
                                   free_blocks,
                                   *volume_._map_index);

  insertBlocksToMap(free_map_, free_blocks);
  init_position_cleared_ = true;
  float res_v = volume_dimension_.cast<float>().x() / volume_resolution_.cast<float>().x();
  // LOG(INFO) << "Planning free_map_  size  " << free_map_.size();

  float step = volume_dimension_.x() / volume_resolution_.x();
  se::exploration::getExplorationPath(discrete_vol_ptr_,
                                      volume_,
                                      free_map_,
                                      frontier_map_,
                                      res_v,
                                      step,
                                      planning_config_,
                                      config_,
                                      pose_ * Tbc_,
                                      path,
                                      cand_views,
                                      exploration_done);
//  std::cout << "[se/denseSLAM] path length " << path.size() <<std::endl;
}

void DenseSLAMSystem::dump_volume(std::string) {

}

void DenseSLAMSystem::renderVolume(unsigned char *out,
                                   const Eigen::Vector2i &outputSize,
                                   int frame,
                                   int raycast_rendering_rate,
                                   const Eigen::Vector4f &k,
                                   float largestep) {

  if (frame % raycast_rendering_rate == 0) {
    const float step = volume_dimension_.x() / volume_resolution_.x();
    renderVolumeKernel(volume_,
                       out,
                       outputSize,
                       *(this->viewPose_) * getInverseCameraMatrix(k),
                       nearPlane,
                       farPlane * 2.0f,
                       mu_,
                       step,
                       largestep,
                       this->viewPose_->topRightCorner<3, 1>(),
                       ambient,
                       !(this->viewPose_->isApprox(raycast_pose_)),
                       vertex_,
                       normal_);
  }
}

void DenseSLAMSystem::renderTrack(unsigned char *out, const Eigen::Vector2i &outputSize) {
  renderTrackKernel(out, tracking_result_.data(), outputSize);
}

void DenseSLAMSystem::renderDepth(unsigned char *out, const Eigen::Vector2i &outputSize) {
  renderDepthKernel(out, float_depth_.data(), outputSize, nearPlane, farPlane);
}

void DenseSLAMSystem::dump_mesh(const std::string filename) {

  std::vector<Triangle> mesh;
  auto inside = [](const Volume<FieldType>::value_type &val) {
    // meshing::status code;
    // if(val.y == 0.f)
    //   code = meshing::status::UNKNOWN;
    // else
    //   code = val.x < 0.f ? meshing::status::INSIDE : meshing::status::OUTSIDE;
    // return code;
    // std::cerr << val.x << " ";
    return val.x < 0.f;
  };

  auto select = [](const Volume<FieldType>::value_type &val) {
    return val.x;
  };

  se::algorithms::marching_cube(*volume_._map_index, select, inside, mesh);
  writeVtkMesh(filename.c_str(), mesh);
}

bool DenseSLAMSystem::getFrontierVoxelMap(map3i &frontier_map) {
  frontier_map = frontier_map_;
  return true;
}
