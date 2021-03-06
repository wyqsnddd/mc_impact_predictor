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
  int nOsdEe = getOsd_()->getEeNum();
  int dim = getOsd_()->getJacobianDim();

  // Arnaud: this assert fails
  //std::cout<<"nOsdEe: "<< nOsdEe << ", numEe_: "<<numEe_<<std::endl;
  //assert(nOsdEe <= numEe_);

  A_.resize(nOsdEe * dim, dof + dim * (nOsdEe));
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

  // Arnaud this assignement has the wrong size!
  // (9x52) <- (9x9)
  A_.block(0, dof, dim * nOsdEe, dim * nOsdEe) = -getOsd_()->getLambdaMatrixInv();
}

} // namespace mc_impact
