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
#include <iostream>
#include <map>

namespace mc_impact
{

struct qpEstimatorParameter
{
  double Qweight = 20;
  // std::vector<std::string> impactBodyNames={"r_wrist"};
  std::map<std::string, Eigen::Vector3d> impactNameAndNormals;
  double impactDuration = 0.005;
  double timeStep = 0.005;
  double coeFrictionDeduction = 0.2;
  double coeRes = 0.8;
  int dim = 3;
  bool useLagrangeMultiplier = false;
  bool useJsd = true;
  bool useOsd = true;
};

struct endEffector
{
  // size_t startingNumber; // Starting number in the vector
  int uniqueIndex;
  Eigen::Vector3d eeVJump;
  Eigen::Vector3d estimatedImpulse;
  Eigen::Vector3d estimatedAverageImpulsiveForce;
  Eigen::Vector3d checkForce;
  Eigen::MatrixXd jacobianDeltaF;

  // The other impulses transformed to the local contact frame
  sva::ForceVecd perturbedWrench;
};

} // namespace mc_impact
