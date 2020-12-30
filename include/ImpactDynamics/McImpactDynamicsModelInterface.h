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
  
    std::cout<<RobotInterface::info<<"Destructing McTwoDimModelBridge." <<RobotInterface::reset<<std::endl;
  }

  inline void setHostCtl(mc_control::fsm::Controller * ctlPtr)
  {

    if(hostCtlPtr_ == nullptr)
    {
      hostCtlPtr_ = ctlPtr;
    }
    else
    {
      throw std::runtime_error("The host fsm controller is already set!");
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
      throw std::runtime_error("The host fsm controller is not set!");
    }
  }
	
 
}; // end of the McTwoDimModelBridge

} // namespace mc_impact
