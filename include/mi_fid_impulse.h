#pragma once

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotModule.h>

#include <RBDyn/Momentum.h>


#include "ImpactDynamicsModelInterface.h"
#include "mi_equality.h"
#include "mi_impactModel.h"
#include "mi_osd.h"
#include "mi_utils.h"
#include <assert.h>
#include <map>

namespace mc_impact
{

class mi_fid_impulse: public mi_equality

/** \brief Equalize the centroidal momentum jump to the sum of the predicted impulse.
 */

{
public:
  /*!
      \param robot reference to the robot model used by the QP
      \param linearJacobian whether use the linear part of the Jacobian?
      \param psdPtr pointer to the OSD model, where we extract the contact information
      */
  mi_fid_impulse(const std::shared_ptr<mi_osd> osdPtr,
		  const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels,
		  const std::shared_ptr<rbd::CentroidalMomentumMatrix> cmmPtr,
		  const std::shared_ptr<mc_impact::TwoDimModelBridge> twoDimmodel 
  );

  ~mi_fid_impulse() {}

  inline std::string nameEq() const override
  {
    return "FIDImpulseConstraint";
  }

  /*!
   * This needs to be called in every iteration only once
   */
  void update() override;

  inline std::shared_ptr<rbd::CentroidalMomentumMatrix> getCmm() const
  {
    return cmmPtr_;
  }

  inline std::shared_ptr<mc_impact::TwoDimModelBridge> getFidModel() const
  {
    return twoDimFidModelPtr_;
  }

private:
  void reset_() override;

  const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels_;

  std::shared_ptr<rbd::CentroidalMomentumMatrix> cmmPtr_;

  std::shared_ptr<mc_impact::TwoDimModelBridge> twoDimFidModelPtr_;
};
} // namespace mc_impact
