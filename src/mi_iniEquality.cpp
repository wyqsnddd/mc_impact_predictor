# include "mi_iniEquality.h"


mi_iniEquality::mi_iniEquality(
		const std::shared_ptr<mi_osd> & osdPtr,
		const std::shared_ptr<mi_impactModel> & impactPtr,
		const int & numEe
		) :  mi_equality(osdPtr), impactPtr_(impactPtr), numEe_(numEe)
{
  reset_();
}

void mi_iniEquality::reset_()
{
  int dof = getOsd_()->getDof();
  //int nContactEe = static_cast<int>(getOsd_()->getContactNum());
  int dim = getOsd_()->getJacobianDim(); 

  A_.resize(dim, dof + dim*(numEe_) );
  A_.setZero();
  b_.resize(dim);
  b_.setZero();
}

void mi_iniEquality::update()
{
  int dof = getOsd_()->getDof();
  int dim = getOsd_()->getJacobianDim(); 
  A_.block(0, 0, dim, dof) = getImpactModel()->getJacobian();
  b_ = getImpactModel()->getEeVelocityJump();

}

