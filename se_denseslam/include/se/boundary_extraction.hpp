/**
 * Information-theoretic exploration, OMPL Motion Validator Skeleton.
 *
 * Copyright (C) 2019 Imperial College London.
 * Copyright (C) 2019 ETH Zürich.
 *
 * @file boundary_extraction.hpp
 * @author Anna Dai
 * @date August 22, 2019
 */

#ifndef SUPEREIGHT_BOUNDARY_EXTRACTION_HPP
#define SUPEREIGHT_BOUNDARY_EXTRACTION_HPP

#include <cstdlib>
#include <map>

#include "continuous/volume_template.hpp"
#include "se/octree.hpp"
#include "se/functors/data_handler.hpp"
#include "se/node_iterator.hpp"
#include "se/utils/eigen_utils.h"

template<typename T> using Volume = VolumeTemplate<T, se::Octree>;

static inline void insertBlocksToSet(set3i &blocks_map, set3i *blocks) {
  if (blocks->size() == 0) return;
  for (auto it = blocks->begin(); it != blocks->end(); ++it) {
    blocks_map.emplace(*it);
  }
}

static inline void removeFromSet(set3i &blocks_map, set3i *blocks_to_remove){
  if(blocks_to_remove->size() == 0 || blocks_map.size()==0) return;
  for(const auto & val : *blocks_to_remove){
     blocks_map.erase(val);
  }
}

/**
* used for sphere initialization
*/
static inline void insertBlocksToMap(set3i &blocks_map, mapvec3i *blocks) {
  if (blocks->size() == 0) return;
  for (auto it = blocks->begin(); it != blocks->end(); ++it) {
    blocks_map.emplace(it->first);
  }
}

// level at leaf level
/**
 * @brief get the boundaries from the free map
 * @param octree_ptr_
 * @param blocks_map   morton code of voxelblocks with free voxels
 * @param lower_bound[out]
 * @param upper_bound[out]
 */
template<typename T>
static inline void getFreeMapBounds(const std::shared_ptr<se::Octree<T> > octree_ptr_,
                                    const set3i &blocks_map,
                                    Eigen::Vector3i &lower_bound,
                                    Eigen::Vector3i &upper_bound) {


  auto it_beg = blocks_map.begin();
  auto it_end = blocks_map.end();
  Eigen::Vector3i lower_bound_tmp;
  Eigen::Vector3i upper_bound_tmp;
  se::Node<T> *node = nullptr;
  bool is_block = false;


  se::key_t lower_bound_morton = *it_beg;
  Eigen::Vector3i lower_block_coord = se::keyops::decode(lower_bound_morton);
  octree_ptr_->fetch_octant(lower_block_coord(0), lower_block_coord(1), lower_block_coord(2), node, is_block);
  bool valid_lower = is_block;
  --it_end;
  se::key_t upper_bound_morton = *it_end;
  Eigen::Vector3i upper_block_coord = se::keyops::decode(upper_bound_morton);
  octree_ptr_->fetch_octant(upper_block_coord(0), upper_block_coord(1), upper_block_coord(2), node, is_block);
  bool valid_upper = is_block;

  // std::cout << "upper " << upper_block_coord.format(InLine) << " lower "<< lower_block_coord.format(InLine)
  // << std::endl;
  lower_bound = upper_block_coord;
  upper_bound = lower_block_coord;
  while (it_beg != it_end) {
    ++it_beg;
    lower_bound_morton = *it_beg;
    lower_bound_tmp = se::keyops::decode(lower_bound_morton);
    octree_ptr_->fetch_octant(lower_bound_tmp(0), lower_bound_tmp(1), lower_bound_tmp(2), node, is_block);
    valid_lower = is_block;

    --it_end;
    upper_bound_morton = *it_end;
    upper_bound_tmp = se::keyops::decode(upper_bound_morton);
    // upper_bound_tmp += Eigen::Vector3i(8,8,8);
    octree_ptr_->fetch_octant(upper_bound_tmp(0), upper_bound_tmp(1), upper_bound_tmp(2), node, is_block);
    valid_upper = is_block;
    // std::cout << "upper " << upper_bound_tmp.format(InLine) << " lower "<< lower_bound_tmp.format(InLine)
    // << std::endl;
    if(valid_lower){
      lower_bound.x() = lower_bound_tmp.x() < lower_bound.x() ? lower_bound_tmp.x() : lower_bound.x();
      lower_bound.y() = lower_bound_tmp.y() < lower_bound.y() ? lower_bound_tmp.y() : lower_bound.y();
      lower_bound.z() = lower_bound_tmp.z() < lower_bound.z() ? lower_bound_tmp.z() : lower_bound.z();
    }
    if (valid_upper){
      upper_bound.x() = upper_bound_tmp.x() +1> upper_bound.x() ? upper_bound_tmp.x()+1 : upper_bound.x();
      upper_bound.y() = upper_bound_tmp.y() +1> upper_bound.y() ? upper_bound_tmp.y()+1 : upper_bound.y();
      upper_bound.z() = upper_bound_tmp.z() +1> upper_bound.z() ? upper_bound_tmp.z()+1 : upper_bound.z();

    }
  }

}
#endif //SUPEREIGHT_BOUNDARY_EXTRACTION_HPP
