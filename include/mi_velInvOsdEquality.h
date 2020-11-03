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

#include "mi_qpEstimator.h"
#include "mi_equality.h"
#include "mi_impactModel.h"
#include "mi_osd.h"

namespace mc_impact
{

class mi_velInvOsdEquality : public mi_equality
{
public:
  mi_velInvOsdEquality(const std::shared_ptr<mi_osd> osdPtr, 
		  const std::shared_ptr<mi_qpEstimator> qpEstimator,
		  const int & numEe);
  ~mi_velInvOsdEquality() {}

  inline std::string nameEq() const override
  {
    return "InverseOperationSpaceDynamicsImpulseEqualityConstraint";
  }
  void update() override;
  const std::shared_ptr<mi_qpEstimator> getQpEstimator()
  {
    return qpEstimator_;
  }
protected:
  void reset_() override;

  std::shared_ptr<mi_qpEstimator> qpEstimator_;
  
  const int numEe_;
};
} // namespace mc_impact
