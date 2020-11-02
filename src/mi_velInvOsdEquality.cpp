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

#include "mi_velInvOsdEquality.h"

namespace mc_impact
{

mi_velInvOsdEquality::mi_velInvOsdEquality(const std::shared_ptr<mi_osd> & osdPtr,  
		const std::shared_ptr<mi_qpEstimator> qpEstimator,
		const int & numEe)
: mi_equality(osdPtr), qpEstimator_(qpEstimator), numEe_(numEe)
{
  reset_();
}

void mi_velInvOsdEquality::reset_()
{
  int dof = getOsd_()->getDof();
  int nOsdEe = getOsd_()->getEeNum();
  int dim = getOsd_()->getJacobianDim();

  // Arnaud: this assert fails
  // std::cout<<"nOsdEe: "<< nOsdEe << ", numEe_: "<<numEe_<<std::endl;
  // assert(nOsdEe <= numEe_);

  A_.resize(nOsdEe * dim, dof);
  A_.setZero();

  b_.resize(nOsdEe * dim);
  b_.setZero();
}

void mi_velInvOsdEquality::update()
{
  int dof = getOsd_()->getDof();
  //int nOsdEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim();
  /*
    Eigen::VectorXd tempId;
    tempId.resize(dim);
    tempId.setIdentity();
  */
  // Go through all the end-effectors and fill in the Jacobians. 
  for(auto idx = getOsd_()->getEes().begin(); idx != getOsd_()->getEes().end(); ++idx)
  {
    int eeIndex = getOsd_()->nameToIndex_(*idx);
    int location = dim * eeIndex;
    // This is the constraint:
    A_.block(location, 0, dim, dof) = getOsd_()->getJacobian(*idx);

    b_.segment(location, dim) = getQpEstimator()->getImpulse(*idx); 

  }

   b_ = -getOsd_()->getLambdaMatrixInv() * b_;
}

} // namespace mc_impact
