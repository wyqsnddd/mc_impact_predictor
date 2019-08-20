# include "mi_iniEquality.h"


mi_iniEquality::mi_iniEquality(
		const std::shared_ptr<mi_osd> & osdPtr,
		const std::shared_ptr<mi_impactModel> & impactPtr
		) :  mi_equality(osdPtr), impactPtr_(impactPtr)
{
  reset_();
}

void mi_iniEquality::reset_()
{
  int dof = getOsd_()->getDof();
  int nEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim(); 

  A_.resize(dim, dof + nEe*dim);
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

