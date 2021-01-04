#pragma once

#include <mc_control/fsm/Controller.h>
#include <RBDyn/Momentum.h>

#include "ImpactDynamicsModelInterface.h"


namespace mc_impact{


class McTwoDimModelBridge : public TwoDimModelBridge 
/*! \brief we add mc_rtc compatible log-entries and GUI components
 */
{
 public: 
  McTwoDimModelBridge(
		  const std::shared_ptr<RobotInterface::Robot> robotPtr,
                  const ImpactModelParams & params,
		  const TwoDimModelBridgeParams & brigeParams);

  virtual ~McTwoDimModelBridge() {
  
    std::cout<<RoboticsUtils::info<<"Destructing McTwoDimModelBridge." <<RoboticsUtils::reset<<std::endl;
  }

  inline void setHostCtl(mc_control::fsm::Controller * ctlPtr)
  {

    if(hostCtlPtr_ == nullptr)
    {
      hostCtlPtr_ = ctlPtr;
    }
    else
    {
      RoboticsUtils::throw_runtime_error("The host fsm controller is already set!", __FILE__, __LINE__);
    }
  }
  void logImpulseEstimations();

  void removeImpulseEstimations();

 protected:
  mc_control::fsm::Controller * hostCtlPtr_ = nullptr;
  inline mc_control::fsm::Controller * getHostCtl_()
  {
    if(hostCtlPtr_ != nullptr)
    {
      return hostCtlPtr_;
    }
    else
    {
      RoboticsUtils::throw_runtime_error("The host fsm controller is not set!", __FILE__, __LINE__);
    }
  }
	
 
}; // end of the McTwoDimModelBridge

} // namespace mc_impact
