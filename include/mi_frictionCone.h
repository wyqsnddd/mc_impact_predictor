#pragma once

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotModule.h>

#include <RBDyn/Momentum.h>


#include "mi_inEquality.h"
#include "mi_osd.h"
#include "mi_utils.h"
#include <assert.h>
#include <map>

namespace mc_impact
{

class mi_frictionCone: public mi_inEquality 

/** \brief Restrict the impulse inside the friction cone.  
 */

{
public:
  /*!
      \param robot reference to the robot model used by the QP
      \param linearJacobian whether use the linear part of the Jacobian?
      \param psdPtr pointer to the OSD model, where we extract the contact information
      */
  mi_frictionCone(const std::shared_ptr<mi_osd> osdPtr);

  ~mi_frictionCone() {}

  inline std::string nameIeq() const override
  {
    return "FrictionConeConstraint";
  }

  /*!
   * This needs to be called in every iteration only once
   */
  void update() override;

private:
  void reset_() override;

  double getMiu_() const
  {
    return miu_; 
  }
  double miu_; ///< The friction coefficient

};
} // namespace mc_impact
