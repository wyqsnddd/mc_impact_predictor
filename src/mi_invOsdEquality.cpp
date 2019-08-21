#include "mi_invOsdEquality.h"

namespace mc_impact
{

mi_invOsdEquality::mi_invOsdEquality(const std::shared_ptr<mi_osd> & osdPtr, const int & numEe)
: mi_equality(osdPtr), numEe_(numEe)
{
  reset_();
}

void mi_invOsdEquality::reset_()
{
  int dof = getOsd_()->getDof();
  int nOsdEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim();

  assert(nOsdEe <= numEe_);

  A_.resize(nOsdEe * dim, dof + dim * (numEe_));
  A_.setZero();

  // b_ will stay zero
  b_.resize(nOsdEe * dim);
  b_.setZero();
}

void mi_invOsdEquality::update()
{
  int dof = getOsd_()->getDof();
  int nOsdEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim();
  /*
    Eigen::VectorXd tempId;
    tempId.resize(dim);
    tempId.setIdentity();
  */
  // Go through all the end-effectors
  for(auto idx = getOsd_()->getEes().begin(); idx != getOsd_()->getEes().end(); ++idx)
  {
    int eeIndex = getOsd_()->nameToIndex_(*idx);
    int location = dim * eeIndex;
    // This is the constraint:
    A_.block(location, 0, dim, dof) = getOsd_()->getJacobian(*idx);

    // A_.block(location, 0, dim, dof) = getOsd_()->getEffectiveLambdaMatrix(*idx);
    // A_.block(location, dof + location, dim, dim) = -tempId;
  }

  A_.block(0, dof, dim * nOsdEe, dim * nOsdEe) = -getOsd_()->getLambdaMatrixInv();
}

} // namespace mi_impact
