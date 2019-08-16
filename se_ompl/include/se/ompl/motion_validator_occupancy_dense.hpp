/**
 * Probabilistic Trajectory Planning, OMPL Motion Validator Dense.
 *
 * Copyright (C) 2018 Imperial College London.
 * Copyright (C) 2018 ETH Zürich.
 *
 * @file MotionValidatorOccupancyDense.hpp
 * @author Nils Funk
 * @date September 14, 2017
 */

#ifndef EXPLORATION_MOTIONVALIDATOROCCUPANCYDENSE_HPP
#define EXPLORATION_MOTIONVALIDATOROCCUPANCYDENSE_HPP

/**
 * Motion Planning, OMPL Motion Validator.
 *
 * Copyright (C) 2018 Imperial College London.
 * Copyright (C) 2018 ETH Zürich.
 *
 * @file MotionValidator.hpp
 *
 * @date June 03, 2017
 */

#include <glog/logging.h>
#include <ctime>
#include <iostream>
#include <stdexcept>

#include "ompl/base/MotionValidator.h"
#include "ompl/base/SpaceInformation.h"

#include "prob_collision_checker.hpp"
#include "se/occupancy_world.hpp"
#include "se/utils/ompl_to_eigen.hpp"

namespace ob = ompl::base;

namespace se {
namespace exploration {

/**
 * Class for interfacing the collsion checking with OMPL.
 *
 * @note Derived from class ompl::base::MotionValidator.
 * @warning
 */

template<typename FieldType>
class MotionValidatorOccupancyDense : public ompl::base::MotionValidator {
 public:
  MotionValidatorOccupancyDense(const ompl::base::SpaceInformationPtr &si,
                                const std::shared_ptr<ProbCollisionChecker<FieldType> > &pcc,
                                const double min_flight_corridor_radius)
      :
      ompl::base::MotionValidator(si),
      pcc_(pcc),
      min_flight_corridor_radius_(min_flight_corridor_radius) {

    stateSpace_ = si_->getStateSpace().get();
    if (stateSpace_ == nullptr)
      throw std::runtime_error("No state space for motion validator");
  }

  /**
   * Checks if the current robot state is valid.
   * @param [in] state The current robot state.
   * @return True if valid, false otherwise.
   */
  virtual bool checkMotion(const ompl::base::State *s1,
                           const ompl::base::State *s2) const override {
    if (!si_->satisfiesBounds(s2)) {
      invalid_++;
      return false;
    }

    Eigen::Vector3d start = OmplToEigen::convertState(*s1);
    Eigen::Vector3d ending = OmplToEigen::convertState(*s2);

    if (pcc_->checkSegmentFlightCorridor(start, ending, 0, min_flight_corridor_radius_)) {
      return true;
    }

    invalid_++;
    return false;
  }

  bool checkMotion(const ompl::base::State *s1,
                   const ompl::base::State *s2,
                   std::pair<ompl::base::State *, double> &lastValid) const override {
    if (!si_->satisfiesBounds(s2)) {
      invalid_++;
      return false;
    }

    /* assume motion starts in a valid configuration so s1 is valid */
    int nd = stateSpace_->validSegmentCount(s1, s2);

    /* temporary storage for the checked state */
    ob::State *test = si_->allocState();
    ob::State *test_prev = si_->allocState();

    for (int j = 1; j <= nd; ++j) {
      stateSpace_->interpolate(s1, s2, (double) j / (double) nd, test);
      stateSpace_->interpolate(s1, s2, (double) (j - 1) / (double) nd, test_prev);

      Eigen::Vector3d start = OmplToEigen::convertState(*test_prev);
      Eigen::Vector3d ending = OmplToEigen::convertState(*test);

      if (!pcc_->checkSegmentFlightCorridor(start, ending, 0, min_flight_corridor_radius_)) {
        lastValid.second = (double) (j - 1) / (double) nd;
        if (lastValid.first != nullptr)
          stateSpace_->interpolate(s1, s2, lastValid.second, lastValid.first);
        invalid_++;
        si_->freeState(test);
        si_->freeState(test_prev);
        return false;
      }
      si_->freeState(test);
      si_->freeState(test_prev);
    }

    valid_++;
    return true;
  }

 private:
  std::shared_ptr<ProbCollisionChecker<FieldType> > pcc_ = nullptr;
  double min_flight_corridor_radius_;
  ob::StateSpace *stateSpace_;
};

} // namespace exploration
}
#endif //EXPLORATION_MOTIONVALIDATOROCCUPANCYDENSE_HPP