/* Copyright 2019 CNRS-UM LIRMM
 *
 * \author Yuquan Wang, Arnaud Tanguy
 *
 *
 *
 * mc_impact_predictor is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * mc_impact_predictor is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mc_impact_predictor. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <mc_rbdyn/Robots.h>

#include <Eigen/Dense>
#include <assert.h>
#include <iostream>
#include <limits>
#include <map>

namespace mc_impact
{

#ifndef COLOUR_PRINT
#  define COLOUR_PRINT
const std::string red("\033[0;31m");
const std::string green("\033[1;32m");
const std::string yellow("\033[1;33m");
const std::string cyan("\033[0;36m");
const std::string magenta("\033[0;35m");
const std::string reset("\033[0m");
#endif

struct qpEstimatorParameter
{
  double Qweight = 20;
  // std::vector<std::string> impactBodyNames={"r_wrist"};
  std::map<std::string, Eigen::Vector3d> impactNameAndNormals;
  double impactDuration = 0.005;
  double qWeightLambda = 0.005;
  double timeStep = 0.005;
  double coeFrictionDeduction = 0.2;
  double coeRes = 0.8;
  int dim = 3;
  bool useLagrangeMultiplier = false;
  bool impactModelBodyJacobian = true;
  bool osdBodyJacobian = true;
  bool useJsd = false;
  bool useOsd = false;
  bool useImpulseBalance = false;
  bool useFidImpulse = false;
  bool useContactConstraint = false;
  // The inequalities: 
  bool useUnilateralContactConstraint = false;
  bool useFrictionCone = false; ///< Restrict the impulse inside the friction cone.

  bool testWeightedQp = false;
  // 0: The default option which minimize the sum of momentum
  // 1: Momentum conservation using the Spatial Jacobian
  // 2: Momentum conservation using the Body Jacobian
  int objectiveChoice = 0;
  std::string name = "qp Estimator";
};

struct endEffector
{
  // size_t startingNumber; // Starting number in the vector
  int uniqueIndex;
  Eigen::Vector3d eeVJump;
  Eigen::Vector3d estimatedImpulse;
  Eigen::Vector3d estimatedAverageImpulsiveForce;
  Eigen::Vector3d checkForce;
  Eigen::Vector3d rssForce;

  Eigen::MatrixXd jacobianDeltaF;

  // The other impulses transformed to the local contact frame
  sva::ForceVecd perturbedWrench;
};

enum class ImpactModelType
{
  OneDimModel,
  TwoDimModel
};
enum class TwoDimModelCase
{
  PushWall
};

struct ImpactModelParams
{
  std::string iBodyName;
  Eigen::Vector3d inertial_surfaceNormal = Eigen::Vector3d::Zero();
  bool useBodyJacobian = true;
  double iDuration = 0.005;
  double timeStep = 0.005;
  double coeF = 0.2;
  double coeR = 0.8;
  int dim = 3;
  ImpactModelType impactModel = ImpactModelType::OneDimModel;
  Eigen::Vector3d eePosition = Eigen::Vector3d::Zero(); ///< Impacting end-effector position 
}; // End of the ImpactModelParams.

inline bool areSame(double a, double b)
{
  return fabs(a - b) < 1e-3;
}

} // namespace mc_impact
