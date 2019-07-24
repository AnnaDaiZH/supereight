/*
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
#ifndef DATA_HANDLER_HPP
#define DATA_HANDLER_HPP
#include "../utils/math_utils.h"
#include "../node.hpp"
#include "../octree.hpp"
#include "../neighbors/neighbor_gather.hpp"

template<typename SpecialisedHandlerT, typename NodeT>
class DataHandlerBase {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typename NodeT::value_type get() {
    return static_cast<SpecialisedHandlerT *>(this)->get();
  }
  void set(const typename NodeT::value_type &val) {
    static_cast<SpecialisedHandlerT *>(this)->set(val);
  }

  virtual Eigen::Vector3i getNodeCoordinates() {
  };

  virtual void occupancyUpdated(const bool o) {
  };

  virtual bool occupancyUpdated() {
  };
  virtual bool isFrontier() {
  };

  virtual key_t get_morton_code() {};

  virtual int get_child_idx() {};

  virtual NodeT *get_node() {};

};

template<typename FieldType>
class VoxelBlockHandler : DataHandlerBase<VoxelBlockHandler<FieldType>,
                                          se::VoxelBlock<FieldType> > {

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> vec3i;

  VoxelBlockHandler(se::VoxelBlock<FieldType> *ptr, Eigen::Vector3i v) : _block(ptr), _voxel(v) {}

  typename se::VoxelBlock<FieldType>::value_type get() {
    return _block->data(_voxel);
  }

  void set(const typename se::VoxelBlock<FieldType>::value_type &val) {
    _block->data(_voxel, val);
  }

  Eigen::Vector3i getNodeCoordinates() {
    return _block->coordinates();
  }

  void occupancyUpdated(const bool o) {
    _block->occupancyUpdated(o);
  }

  bool occupancyUpdated() {
    return _block->occupancyUpdated();
  }

  bool isFrontier(const se::Octree<FieldType> &map) {
    vec3i face_neighbour_voxel(6);
    face_neighbour_voxel[0] << _voxel.x() - 1, _voxel.y(), _voxel.z();
    face_neighbour_voxel[1] << _voxel.x() + 1, _voxel.y(), _voxel.z();
    face_neighbour_voxel[2] << _voxel.x(), _voxel.y() - 1, _voxel.z();
    face_neighbour_voxel[3] << _voxel.x(), _voxel.y() + 1, _voxel.z();
    face_neighbour_voxel[4] << _voxel.x(), _voxel.y(), _voxel.z() - 1;
    face_neighbour_voxel[5] << _voxel.x(), _voxel.y(), _voxel.z() + 1;

    for (const auto &face_voxel : face_neighbour_voxel) {

      // check if face voxel is inside same voxel block
      if ((_voxel.x() / BLOCK_SIDE) == (face_voxel.x() / BLOCK_SIDE)
          && (_voxel.y() / BLOCK_SIDE) == (face_voxel.y() / BLOCK_SIDE)
          && (_voxel.z() / BLOCK_SIDE) == (face_voxel.z() / BLOCK_SIDE)) {
        // same voxel block
        if (_block->data(face_voxel).st == voxel_state::kUnknown)
          return true;
      } else {
        // not same voxel block => check if neighbour is a voxel block
        se::Node<FieldType> *node = nullptr;
        bool is_voxel_block;
        map.fetch_octant(face_voxel(0), face_voxel(1), face_voxel(2), node, is_voxel_block);

        if (is_voxel_block) {
          // neighbour is a voxel block
          se::VoxelBlock<FieldType> *block = static_cast<se::VoxelBlock<FieldType> *> (node);
          if (block->data(face_voxel).st == voxel_state::kUnknown)
            return true;
//          return block->data(face_voxel).st == voxel_state::kUnknown;
        } else {
          // neighbour is a node
          if (map.get(se::keyops::decode(node->code_)).st == voxel_state::kUnknown)
            return true;
        }
      }
    }
    return false;
  }

  key_t get_morton_code() {
    return _block->code_;
  }
// only for nodes
  int get_child_idx() {
    return -1;
  };

  se::VoxelBlock<FieldType> *get_node() {
    return _block;
  }

 private:
  se::VoxelBlock<FieldType> *_block;
  Eigen::Vector3i _voxel;
};

template<typename FieldType>
class NodeHandler : DataHandlerBase<NodeHandler<FieldType>, se::Node<FieldType> > {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  NodeHandler(se::Node<FieldType> *ptr, int i) : _node(ptr), _idx(i) {}

  typename se::Node<FieldType>::value_type get() {
    return _node->value_[_idx];
  }

  void set(const typename se::Node<FieldType>::value_type &val) {
    _node->value_[_idx] = val;
  }

  Eigen::Vector3i getNodeCoordinates() {
    return  _node->childCoordinates(_idx);
  }

  void occupancyUpdated(const bool o) {
    _node->occupancyUpdated(o);
  }

  bool occupancyUpdated() {
    return _node->occupancyUpdated();
  }
  /*
   * check if neighbour octant are unknown
   */
  bool isFrontier(const se::Octree<FieldType> &map) {
    // find out at which level the node is
    key_t morton_parent = _node->code_;
    Eigen::Vector3i coord = _node->childCoordinates(_idx);
    int level_parent = se::keyops::level(morton_parent);
    int side_parent = _node->side_;
    int level_leaf = level_parent;
    while (side_parent != BLOCK_SIDE) {
      side_parent= side_parent >> 1;
      level_leaf++;
    }
    int scale = level_leaf- level_parent;
    // when this node is a leaf node , we only want to have frontiers at the higher resolution /
    // voxel level
    // how ot check if the node is a voxel block
    // NILS: if node at leaf level : node = voxel block
    // get face octants
//    std::cout << "curr node status " << map.get(coord).st << " by val " << _node->value_[_idx].st
//    <<"prob "<<   _node->value_[_idx].x <<
//    std::endl;
//    std::cout << " level leaf " << level_leaf << " side parent after " << side_parent << " scale "
//              << scale << std::endl;

    for (size_t i = 0; i < 6; ++i) {
      // Compute the neighbor octant coordinates.
      const int neighbor_x = coord.x() + scale * BLOCK_SIDE * se::face_neighbor_offsets[i].x();
      const int neighbor_y = coord.y() + scale * BLOCK_SIDE * se::face_neighbor_offsets[i].y();
      const int neighbor_z = coord.z() + scale * BLOCK_SIDE * se::face_neighbor_offsets[i].z();
      // neighbour node is allocated and is unknown
      if (map.fetch_octant(neighbor_x, neighbor_y, neighbor_z, level_parent +1) != nullptr
          && map.get(neighbor_x, neighbor_y, neighbor_z).st == voxel_state::kUnknown) {
//        std::cout << "[se/datahandler] " << coord.format(InLine) << " is frontier " << neighbor_x
//                  << " " << neighbor_y << " " << neighbor_z << std::endl;
        return true;
      }
    }
    return false;
  }
  // returns the morton code of this node.

  key_t get_morton_code() {
    const int level = se::keyops::level(_node->code_) + 1; // level of the node which is being
    // updated
    // TODO find a way to pass max_level around or make it constant
    const int max_level = 7;
    const Eigen::Vector3i octant_coord = getNodeCoordinates();
    // compute morton code of current node / voxelblock
    const key_t morton_code = se::keyops::encode(octant_coord.x(),
                                                octant_coord.y(),
                                                octant_coord.z(),
                                                level,
                                                max_level);
    return morton_code;
  };

  int get_child_idx() {
    return _idx;
  }

  se::Node<FieldType> *get_node() {
    return _node;
  }


 private:
  se::Node<FieldType> *_node; // parent node
  int _idx; //  index of this node
};

#endif
