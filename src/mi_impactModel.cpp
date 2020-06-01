#include "mi_impactModel.h"

namespace mc_impact
{

void mi_impactModel::update(const Eigen::Vector3d & surfaceNormal)
{
  local_surfaceNormal_ = surfaceNormal;
  update_();
}

void mi_impactModel::update()
{

  if(useBodyJacobian())
  {
    sva::PTransformd X_0_ee = simRobot_.bodyPosW(getImpactBody());
    local_surfaceNormal_ = X_0_ee.rotation() * inertial_surfaceNormal_;
    local_surfaceNormal_.normalize();
  }
  else{
    local_surfaceNormal_ = inertial_surfaceNormal_;
  }
  update_();
}

void mi_impactModel::updateJacobian_()
{
  Eigen::MatrixXd tempJacobian;

  if(useBodyJacobian()){
    tempJacobian  = jacPtr_->bodyJacobian(simRobot_.mb(), simRobot_.mbc());
  }else{
    tempJacobian  = jacPtr_->jacobian(simRobot_.mb(), simRobot_.mbc());
  }

  jacPtr_->fullJacobian(simRobot_.mb(), tempJacobian.block(3, 0, 3, tempJacobian.cols()), jacobian_);
}

void mi_impactModel::update_()
{

  updateJacobian_();

  // std::cout<<"impactModel: surfaceNormal is: "<<surfaceNormal.transpose()<<std::endl;
  Eigen::Matrix3d tempProjector = local_surfaceNormal_ * local_surfaceNormal_.transpose();
  // Eigen::Matrix3d tempProjector;
  // tempProjector.setIdentity();
  // Eigen::Matrix3d tempNullProjector = Eigen::Matrix3d::Identity() - tempProjector;

  // Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes()) * tempProjector );

  // std::cout<<"impactModel: temp projector is: "<<std::endl<<tempReductionProjector<<std::endl;
  // Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes()) * tempProjector + getCoeFricDe() * tempNullProjector);
  Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes()) * tempProjector);
  // reductionProjector_ = -((1 + getCoeRes()) * tempProjector + getCoeFricDe() *
  // tempNullProjector)*osdPtr_->getJacobian(getImpactBody());

  // reductionProjector_ = tempReductionProjector*osdPtr_->getJacobian(getImpactBody());
  reductionProjector_ = tempReductionProjector * getJacobian();

  //temp_q_vel_ = (alpha + alphaD * getTimeStep());
  robotJointVel_ = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alphaD);
  // eeV_ = osdPtr_->getJacobian(getImpactBody()) *  temp_q_vel_;
  eeV_ = getJacobian() * getJointVel();

  contactVel_ = tempProjector * eeV_;
  /*
    std::cout<<"impactModel: q_vel is: "<<std::endl<<alpha.transpose()<<std::endl;
    std::cout<<"impactModel: q_acc is: "<<std::endl<<alphaD.transpose()<<std::endl;
    std::cout<<"impactModel: time step is: "<<getTimeStep()<<std::endl;

    std::cout<<"impactModel: jacobian test is: "<<std::endl<<osdPtr_->getJacobian("r_wrist")*
    osdPtr_->getJacobian("r_wrist").transpose()<<std::endl; std::cout<<"impactModel: qp velocity is:
    "<<eeV_.transpose()<<std::endl;
    */
  // deltaV_ = reductionProjector_ * temp_q_vel_ ;
  deltaV_ = tempReductionProjector * eeV_;

  // std::cout<<"impactModel: deltaV is: "<<deltaV_.transpose()<<std::endl;
}

} // namespace mc_impact
