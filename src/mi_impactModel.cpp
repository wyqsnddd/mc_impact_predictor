# include "mi_impactModel.h"

void mi_impactModel::update(const Eigen::Vector3d & surfaceNormal)
{
  Eigen::Matrix3d tempProjector = surfaceNormal * surfaceNormal.transpose();
  Eigen::Matrix3d tempNullProjector = Eigen::Matrix3d::Identity() - surfaceNormal * surfaceNormal.transpose();


  Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes()) * tempProjector + getCoeFricDe() * tempNullProjector);
  Eigen::VectorXd alpha = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alphaD);

  deltaV_ = tempReductionProjector * osdPtr_->getJacobian(getImpactBody()) * (alpha + alphaD * getImpactDuration());

}


