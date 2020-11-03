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

#include "mi_velIniEquality.h"

namespace mc_impact
{

mi_velIniEquality::mi_velIniEquality(const std::shared_ptr<mi_osd> osdPtr,
                               const std::shared_ptr<mi_impactModel> & impactPtr,
                               const int numEe)
: mi_equality(osdPtr), impactPtr_(impactPtr), numEe_(numEe)
{
  reset_();
}

void mi_velIniEquality::reset_()
{
  int dof = getOsd_()->getDof();
  // int nContactEe = static_cast<int>(getOsd_()->getContactNum());
  int dim = getOsd_()->getJacobianDim();

  A_.resize(dim, dof);
  A_.setZero();
  b_.resize(dim);
  b_.setZero();
}

void mi_velIniEquality::update()
{
  //int dof = getOsd_()->getDof();
  //int dim = getOsd_()->getJacobianDim();
  A_ = getImpactModel()->getJacobian();
  b_ = getImpactModel()->getEeVelocityJump();
}

} // namespace mc_impact
