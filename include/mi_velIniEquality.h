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
#include "ImpactDynamics/OneDimImpactModel.h"
#include "mi_equality.h"

namespace mc_impact
{

class mi_velIniEquality : public mi_equality
/** \brief initial condition given the predicted impact
 * J_k \delta q_dot = kth end-effector velocity jump.
 */
{
public:
  mi_velIniEquality(const std::shared_ptr<mi_osd> osdPtr,
                    const std::shared_ptr<mi_impactModel> & impactPtr,
                    const int numEe);
  ~mi_velIniEquality() {}

  inline std::string nameEq() const override
  {
    return "InitialConditionEqualityConstraint";
  }
  void update() override;
  inline const std::shared_ptr<mi_impactModel> getImpactModel()
  {
    return impactPtr_;
  }

protected:
  void reset_() override;

  const std::shared_ptr<mi_impactModel> impactPtr_;

  const int numEe_;
};
} // namespace mc_impact
