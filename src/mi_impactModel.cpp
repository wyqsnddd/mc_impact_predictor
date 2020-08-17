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

  if(getParams().useBodyJacobian)
  {
    sva::PTransformd X_0_ee = simRobot_.bodyPosW(getParams().iBodyName);
    local_surfaceNormal_ = X_0_ee.rotation() * getParams().inertial_surfaceNormal;
    local_surfaceNormal_.normalize();
  }
  else{
    local_surfaceNormal_ = getParams().inertial_surfaceNormal;
  }
  update_();
}

void mi_impactModel::updateJacobian_()
{
  Eigen::MatrixXd tempJacobian;
  Eigen::MatrixXd tempJacobianDot;

  if(getParams().useBodyJacobian){
    tempJacobian  = jacPtr_->bodyJacobian(simRobot_.mb(), simRobot_.mbc());
    tempJacobianDot  = jacPtr_->bodyJacobianDot(simRobot_.mb(), simRobot_.mbc());
  }else{
    tempJacobian  = jacPtr_->jacobian(simRobot_.mb(), simRobot_.mbc());
    tempJacobianDot  = jacPtr_->jacobianDot(simRobot_.mb(), simRobot_.mbc());
  }

  jacPtr_->fullJacobian(simRobot_.mb(), tempJacobian.block(3, 0, 3, tempJacobian.cols()), jacobian_);

  jacPtr_->fullJacobian(simRobot_.mb(), tempJacobianDot.block(3, 0, 3, tempJacobian.cols()), jacobianDot_);
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
  
  Eigen::Matrix3d tempReductionProjector = -((1 + getParams().coeR) * tempProjector);
  // reductionProjector_ = -((1 + getCoeRes()) * tempProjector + getCoeFricDe() *
  // tempNullProjector)*osdPtr_->getJacobian(getImpactBody());

  // reductionProjector_ = tempReductionProjector*osdPtr_->getJacobian(getImpactBody());
  reductionProjector_ = tempReductionProjector * getJacobian();
  reductionProjectorTwo_ = reductionProjector_ + tempReductionProjector * getJacobianDot() * getParams().timeStep;


  //Eigen::VectorXd alpha = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alpha);
  //Eigen::VectorXd alphaD = rbd::dofToVector(simRobot_.mb(), simRobot_.mbc().alphaD);
  // temp_q_vel_ = (alpha + alphaD * getTimeStep());
  //temp_q_vel_ = alpha;
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
