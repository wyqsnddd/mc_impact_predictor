# include "mi_osdEquality.h"


mi_osdEquality::mi_osdEquality(
		const std::shared_ptr<mi_osd> & osdPtr
		) :  mi_equality(osdPtr)
{
  reset_();
}

void mi_osdEquality::reset_()
{
  int dof = getOsd_()->getDof();
  int nEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim(); 

  A_.resize(nEe*dim, dof + nEe*dim);
  A_.setZero();

  // b_ will stay zero
  b_.resize(nEe*dim);
  b_.setZero();
}


void mi_osdEquality::update()
{
  int dof = getOsd_()->getDof();
  //int nEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim(); 

  Eigen::VectorXd tempId;
  tempId.resize(dim);
  tempId.setIdentity();

// Go through all the end-effectors
  for(auto idx = getOsd_()->getEes().begin(); idx!=getOsd_()->getEes().end(); ++idx)
  {
  int eeIndex = getOsd_()->nameToIndex_(*idx);
  int location =  dim*eeIndex;
  // This is the constraint: 
  A_.block(location, 0, dim, dof) = getOsd_()->getEffectiveLambdaMatrix(*idx);
  A_.block(location, dof + location, dim, dim) = -tempId; 
  }
}
