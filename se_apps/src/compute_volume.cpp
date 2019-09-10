/**
 * Copyright (C) 2019 Imperial College London.
 * Copyright (C) 2019 ETH Zürich.
 *
 *
 *
 * @file exploration_progress_evaluation.cpp
 * @author Anna Dai
 * @author Sotiris Papatheodorou
 * @date August 24, 2019
 */



#include <cstdio>
#include <iostream>
#include <string>

#include "se/octree.hpp"
#include "se/post_processing.hpp"



struct arguments {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  std::string filename;
  Eigen::Vector3f map_dim;
  Eigen::Vector3f map_offset;
};



struct arguments parse_arguments(int argc, char** argv) {
  struct arguments args;

  switch (argc) {
    // Map filename and dimensions supplied.
    case 5:
      args.filename = std::string(argv[1]);
      args.map_dim.x() = atof(argv[2]);
      args.map_dim.y() = atof(argv[3]);
      args.map_dim.z() = atof(argv[4]);
	  args.map_offset = Eigen::Vector3f::Constant(0.f);
      break;

    // Map filename, dimensions and offset supplied.
    case 8:
      args.filename = std::string(argv[1]);
      args.map_dim.x()    = atof(argv[2]);
      args.map_dim.y()    = atof(argv[3]);
      args.map_dim.z()    = atof(argv[4]);
      args.map_offset.x() = atof(argv[5]);
      args.map_offset.y() = atof(argv[6]);
      args.map_offset.z() = atof(argv[7]);
      break;

    // Show usage message.
    default:
      std::cout << "Usage: " << argv[0]
          << " FILENAME DIM_X DIM_Y DIM_Z [OFFSET_X OFFSET_Y OFFSET_Z]\n";
      exit(EXIT_FAILURE);
  }

  // Test if the file exists.
  FILE* fp = fopen(args.filename.c_str(), "r");
  if (fp == nullptr) {
    std::cout << "Error: file " << args.filename
        << " does not exist or is not accessible\n";
    exit(EXIT_FAILURE);
  } else {
    fclose(fp);
  }

  return args;
}



int main(int argc, char** argv) {
  // Parse the input arguments.
  const struct arguments args = parse_arguments(argc, argv);

  // Initialize the Octree and load the saved map.
  std::shared_ptr<se::Octree<OFusion> > octree_
      = std::make_shared<se::Octree<OFusion> >();
  octree_->load(args.filename);

  // Get and show Octree info.
  const int   octree_size   = octree_->size();
  const float octree_dim    = octree_->dim();
  const float octree_volume = std::pow(octree_dim, 3.f);
  const float voxel_dim     = octree_->voxelDim();
  const float voxel_volume  = std::pow(voxel_dim, 3.f);
  std::cout << "Octree info -----------------------\n";
  std::cout << "Octree size:   " << octree_size   << " x "
                                 << octree_size   << " x "
                                 << octree_size   << " voxels\n";
  std::cout << "Octree dim:    " << octree_dim    << " x "
                                 << octree_dim    << " x "
                                 << octree_dim    << " m\n";
  std::cout << "Octree volume: " << octree_volume << " m^3\n";
  std::cout << "Voxel dim:     " << voxel_dim     << " m\n";
  std::cout << "Voxel volume:  " << voxel_volume  << " m^3\n";

  // Show map info.
  std::cout << "Map info --------------------------\n";
  std::cout << args.map_dim.x() << " x "
            << args.map_dim.y() << " x "
            << args.map_dim.z() << " m\n";



  // Count the explored voxels.
  float  free_voxel_volume;
  float  occupied_voxel_volume;
  float  free_node_volume;
  float  occupied_node_volume;
  compute_bounded_volume(*octree_, args.map_dim, args.map_offset,
      free_voxel_volume, occupied_voxel_volume,
      free_node_volume, occupied_node_volume);



  // Print results.
  std::cout << "Voxels ----------------------------\n";
  std::cout << "Explored voxel volume: "
      << free_voxel_volume + occupied_voxel_volume << " m^3\n";
  std::cout << "Nodes -----------------------------\n";
  std::cout << "Explored node volume:  "
      << free_node_volume + occupied_node_volume << " m^3\n";
  std::cout << "Results ---------------------------\n";
  std::cout << "Explored volume:       "
      << free_node_volume + occupied_node_volume
      + free_voxel_volume + occupied_voxel_volume << " m^3\n";

  exit(EXIT_SUCCESS);
}

