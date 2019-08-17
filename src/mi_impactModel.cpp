# include "mi_impactModel.h"

void mi_impactModel::update(const Eigen::Vector3d & surfaceNormal)
{

  surfaceNormal_ = surfaceNormal;
  //surfaceNormal_.setOnes(); 

  //std::cout<<"impactModel: surfaceNormal is: "<<surfaceNormal.transpose()<<std::endl;
  Eigen::Matrix3d tempProjector = surfaceNormal_ * surfaceNormal_.transpose();
  //Eigen::Matrix3d tempProjector;
  //tempProjector.setIdentity();
  Eigen::Matrix3d tempNullProjector = Eigen::Matrix3d::Identity() - tempProjector;

  //Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes()) * tempProjector );

  //std::cout<<"impactModel: temp projector is: "<<std::endl<<tempReductionProjector<<std::endl;
  Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes()) * tempProjector + getCoeFricDe() * tempNullProjector);
  //reductionProjector_ = -((1 + getCoeRes()) * tempProjector + getCoeFricDe() * tempNullProjector)*osdPtr_->getJacobian(getImpactBody());
  reductionProjector_ = tempReductionProjector*osdPtr_->getJacobian(getImpactBody());

  Eigen::VectorXd alpha = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alphaD);
  temp_q_vel_ = (alpha + alphaD * getTimeStep());
  eeV_ = osdPtr_->getJacobian(getImpactBody()) *  temp_q_vel_;
/*
  std::cout<<"impactModel: q_vel is: "<<std::endl<<alpha.transpose()<<std::endl;
  std::cout<<"impactModel: q_acc is: "<<std::endl<<alphaD.transpose()<<std::endl;
  std::cout<<"impactModel: time step is: "<<getTimeStep()<<std::endl;

  std::cout<<"impactModel: jacobian test is: "<<std::endl<<osdPtr_->getJacobian("r_wrist")* osdPtr_->getJacobian("r_wrist").transpose()<<std::endl;
  std::cout<<"impactModel: qp velocity is: "<<eeV_.transpose()<<std::endl;
  */
  //deltaV_ = reductionProjector_ * temp_q_vel_ ;
  deltaV_ = tempReductionProjector * eeV_; 

  //std::cout<<"impactModel: deltaV is: "<<deltaV_.transpose()<<std::endl;
}


