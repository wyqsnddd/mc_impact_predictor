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

#include "mi_jsdEquality.h"

namespace mc_impact
{

mi_jsdEquality::mi_jsdEquality(const std::shared_ptr<mi_osd> osdPtr,
                               const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels,
                               const std::map<std::string, endEffector> & endEffectors)
: mi_equality(osdPtr), impactModels_(impactModels), endEffectors_(endEffectors)
{
  /*
  for (auto idx = impactModels_.begin(); idx!= impactModels_.end();++idx)
  {
   impactBodyNames_.push_back(idx->first);
  }
  */
  reset_();

  std::cout << red << "Initialized the joint space dynamics equality constraint" << reset << std::endl;
}

void mi_jsdEquality::reset_()
{
  // Use the endEffectors with contacts
  int dof = getOsd_()->getDof();
  // int nContactEe = static_cast<int>(getOsd_()->getEeNum());
  // int nImpactEe = static_cast<int>(impactModels_.size());
  int nEe = static_cast<int>(endEffectors_.size());
  int dim = getOsd_()->getJacobianDim();

  A_.resize(dof, dof + dim * (nEe));
  A_.setZero();

  // b_ will stay zero
  b_.resize(dof);
  b_.setZero();
}

int mi_jsdEquality::nameToIndex_(const std::string & eeName)
{
  auto tempEe = endEffectors_.find(eeName);
  if(tempEe != endEffectors_.end())
    return tempEe->second.uniqueIndex;
  else
  {
    std::string error_msg =
        std::string("mi_jsdEquality::nameToIndex_: ee-") + eeName + std::string(": does not exist.");
    throw std::runtime_error(error_msg);
  }
}

void mi_jsdEquality::update()
{

  int nRow = getOsd_()->getDof();
  int dim = getOsd_()->getJacobianDim();
  // (0) Fill the mass matrix
  A_.block(0, 0, nRow, nRow) = getOsd_()->getMassMatrix();

  // (1) Fill the contact bodies, which are right after the joint velocity jump

  for(auto idx = getOsd_()->getContactEes().begin(); idx != getOsd_()->getContactEes().end(); ++idx)
  {
    int eeIndex = getOsd_()->nameToIndex_(*idx);
    // int eeIndex = nameToIndex_(*idx);
    A_.block(0, nRow + eeIndex * dim, nRow, dim) = -getOsd_()->getJacobian(*idx).transpose();
  }

  // int osdVarNum = dim*getOsd_()->getEeNum();

  // (2) Fill the impact bodies after the OSD endeffectors

  for(auto idx = impactModels_.begin(); idx != impactModels_.end(); ++idx)
  {
    int eeIndex = nameToIndex_(idx->first);
    // A_.block(0, nRow + osdVarNum + eeIndex*dim, nRow, dim) = - idx->second->getJacobian().transpose();
    A_.block(0, nRow + dim * eeIndex, nRow, dim) = -idx->second->getJacobian().transpose();
  }
}

} // namespace mc_impact
