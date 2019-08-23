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
  // sva::ForceVecd perturbedWrench;
};

} // namespace mc_impact
