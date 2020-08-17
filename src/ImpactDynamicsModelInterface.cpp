#include "ImpactDynamicsModelInterface.h"
namespace mc_impact 
{

ImpactDynamicsModel::ImpactDynamicsModel(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params) : simRobot_(simRobot), params_(params)
{
  std::cout<<green<<"ImpactDynamicsModel is created."<<reset<<std::endl;

}

TwoDimModelBridge::TwoDimModelBridge(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params): ImpactDynamicsModel(simRobot, params)
{

  // Initialize the semiaxes calculator
  
  // Initialize the virtual-contact-point calculator
  
  // Initialize the two-dim-model 
  std::cout<<green<<"TwoDimModelBridge is created."<<reset<<std::endl;
}


void TwoDimModelBridge::update()
{

  // Update the ssa model
  
  // Update the vc model
  vcParams_.com = getRobot().com();
  // Update the twoDim model
  
}

} // End of NameSpace
