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

#include "mi_unilateralContactConstraint.h"

namespace mc_impact
{

mi_unilateralContactConstraint::mi_unilateralContactConstraint(const std::shared_ptr<mi_osd> osdPtr) : mi_inEquality(osdPtr)
{
  
  reset_();

  std::cout << red << "Initialized the unilateral contact inequality constraint" << reset << std::endl;
}

void mi_unilateralContactConstraint::reset_()
{
  // Use the endEffectors with contacts
  int dof = getOsd_()->getDof();
  int nContact = static_cast<int>(getOsd_()->getContactEes().size());

  //int jacDim = getOsd_()->getJacobianDim();
  int constraintDim = 1; 

  A_.resize(constraintDim * nContact, dof);
  A_.setZero();

  // b_ will stay zero
  b_.resize(constraintDim * nContact);
  b_.setZero();
}

void mi_unilateralContactConstraint::update()
{

  int dof = getOsd_()->getDof();
  int constraintDim = 1;

  // Fill the contact bodies 
  int count = 0;
  for(auto idx = getOsd_()->getContactEes().begin(); idx != getOsd_()->getContactEes().end(); ++idx)
  {
    // We use the x-axis row of the Jacobian.
    //A_.block(count * constraintDim, 0, constraintDim, dof) = getOsd_()->getJacobian(*idx).block(0, 0, constraintDim, dof);
    
    // Suppose z is the contact normal direction, otherwise: n^T * Jac
    A_.block(count * constraintDim, 0, constraintDim, dof) = getOsd_()->getJacobian(*idx).block(2, 0, constraintDim, dof);
    count++;
  }
}

} // namespace mc_impact
