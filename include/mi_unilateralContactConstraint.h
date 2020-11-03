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

#pragma once
#include "mi_inEquality.h"
#include "mi_impactModel.h"
#include "mi_osd.h"
#include "mi_utils.h"

namespace mc_impact
{

class mi_unilateralContactConstraint : public mi_inEquality
/** \brief Specify the unilateral contact constraint: contact velocity <= 0
 *  Suppose the contact normal is n 
 *  (n^T * J * \Delta \dot{q})  <= 0
 */
{
public:
  mi_unilateralContactConstraint(const std::shared_ptr<mi_osd> osdPtr);
  ~mi_unilateralContactConstraint() {}

  inline std::string nameIeq() const override
  {
    return "JointSpaceDynamicsImpulseEqualityConstraint";
  }

  /*!
    \param contactEe We use the end-effectors in contact to specify the jsd impulse equation
    */
  void update() override;

protected:
  void reset_() override;
};

} // namespace mc_impact
