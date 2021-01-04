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

#include "mi_velJsdEquality.h"

namespace mc_impact
{

mi_velJsdEquality::mi_velJsdEquality(const std::shared_ptr<mi_osd> osdPtr,
		const std::shared_ptr<mi_qpEstimator> qpEstimator
                               ):mi_equality(osdPtr), qpEstimator_(qpEstimator) 
{
  reset_();

  std::cout << RoboticsUtils::alarm << "Initialized the joint space dynamics equality constraint" << RoboticsUtils::reset << std::endl;
}

void mi_velJsdEquality::reset_()
{
  // Use the endEffectors with contacts
  int dof = getOsd_()->getDof();
  // int nContactEe = static_cast<int>(getOsd_()->getEeNum());
  // int nImpactEe = static_cast<int>(impactModels_.size());

  A_.resize(dof, dof);
  A_.setZero();

  // b_ will stay zero
  b_.resize(dof);
  b_.setZero();
}

void mi_velJsdEquality::update()
{

  // (0) Fill the mass matrix
  A_ = getOsd_()->getMassMatrix();
  b_.setZero();

  // (1) Fill the contact bodies, which are right after the joint velocity jump

  for(auto idx = getOsd_()->getContactEes().begin(); idx != getOsd_()->getContactEes().end(); ++idx)
  {
    //int eeIndex = getOsd_()->nameToIndex_(*idx);
    // int eeIndex = nameToIndex_(*idx);
    b_ += getOsd_()->getJacobian(*idx).transpose() * getQpEstimator()->getImpulse(*idx);
  }

  // int osdVarNum = dim*getOsd_()->getEeNum();

  // (2) Fill the impact bodies after the OSD endeffectors

  for(auto idx = getQpEstimator()->getImpactModels().begin(); idx != getQpEstimator()->getImpactModels().end(); ++idx)
  {
    //int eeIndex = nameToIndex_(idx->first);
    // A_.block(0, nRow + osdVarNum + eeIndex*dim, nRow, dim) = - idx->second->getJacobian().transpose();
    b_ += idx->second->getJacobian().transpose() * getQpEstimator()->getImpulse(idx->first);
  }
}

} // namespace mc_impact
