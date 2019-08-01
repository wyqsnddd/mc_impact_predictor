# include "mi_iniEquality.h"


mi_iniEquality::mi_iniEquality(
		const std::shared_ptr<mi_osd> & osdPtr,
		const mi_impactModel* impactPtr,
		bool dcJacobian
		) :  mi_equality(osdPtr), impactPtr_(impactPtr), dcJacobian_(dcJacobian)
{
  reset_();
}

void mi_iniEquality::reset_()
{
  int dof = getOsd_()->getDof();
  int nEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim(); 
  if(dcJacobian_){
    A_.resize(dof, dof + nEe*dim);
    A_.setZero();
    b_.resize(dof);
    b_.setZero();

  }else{
    A_.resize(dim, dof + nEe*dim);
    A_.setZero();
    b_.resize(dim);
    b_.setZero();
  }
}

void mi_iniEquality::updateDcJacobian_()
{
  int dof = getOsd_()->getDof();
  int dim = getOsd_()->getJacobianDim(); 

  A_.block(0, 0, dof, dof) = Eigen::MatrixXd::Identity(dof, dof); 

  b_ = getOsd_()->getDcJacobianInv(impactPtr_->getImpactBody())*impactPtr_->getEeVelocityJump();

}

void mi_iniEquality::updateJacobian_()
{
  int dof = getOsd_()->getDof();
  int dim = getOsd_()->getJacobianDim(); 
  A_.block(0, 0, dim, dof) = getOsd_()->getJacobian(impactPtr_->getImpactBody());
  b_ = impactPtr_->getEeVelocityJump();

}
void mi_iniEquality::update()
{
  // This is the constraint: 
  //updateJacobain_();
  if(dcJacobian_){
    updateDcJacobian_();
  }else{
    updateJacobian_();
  }
}
