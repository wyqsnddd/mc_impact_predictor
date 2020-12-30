#pragma once
#include <Eigen/Core>

namespace mc_impact{

enum class ImpactModelType
{
  OneDimModel,
  TwoDimModel
};

struct ImpactModelParams
{
  std::string iBodyName;
  std::string iSurfaceName;
  Eigen::Vector3d inertial_surfaceNormal = Eigen::Vector3d::Zero();
  bool useBodyJacobian = true;
  double iDuration = 0.005;
  double timeStep = 0.005;
  double coeF = 0.8;
  double coeR = 0.2;
  int dim = 3;
  ImpactModelType impactModel = ImpactModelType::OneDimModel;
  Eigen::Vector3d eePosition = Eigen::Vector3d::Zero(); ///< Impacting velocity 
}; // End of the ImpactModelParams.



} // namespace mc_impact

