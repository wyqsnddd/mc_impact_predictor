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
  Eigen::Matrix3d inertiaMatrix = Eigen::Matrix3d::Identity();
  ssaPtr_.reset(new FIDynamics::SolveSemiAxes(getRobot().mass(), inertiaMatrix));
  
  // Initialize the virtual-contact-point calculator
  vcPtr_.reset(new FIDynamics::VirtualContactPoint());

  // Initialize the two-dim-model 
  // Energetic coefficient of restitution
  piParams_.e = 0.8;
  // Coefficient of friction 
  piParams_.miu = 0.7;

  twoDimModelPtr_.reset(new FIDynamics::TwoDimModel(piParams_));
  std::cout<<green<<"TwoDimModelBridge is created."<<reset<<std::endl;
}

void TwoDimModelBridge::updatePlanarImpactParams_()
{
  // How to obtain the planar impact model parameters?  
  /*
  params_.tu = ? 
  params_.nu = ? 
  params_.contactPoint = ? 
  */

}

void TwoDimModelBridge::update(const Eigen::Vector3d & impactNormal)
{

  // (0) Compute the whole-body inertia and average velocity
  Eigen::Matrix6d centroidalInertia; 
  centroidalInertia.setIdentity();

  sva::ForceVecd cm;
  sva::ForceVecd av;
  rbd::computeCentroidalInertiaAndVelocity(getRobot().mb(), getRobot().mbc(), getRobot().com(), centroidalInertia, cm, av);

  // (1) Update the ssa model
  
  // Inertia should be the upper corner? check! 
  std::cout<<alarm<< "The centroidal inertia is: " << std::endl << centroidalInertia<< std::endl;
  ssaPtr_->update(getRobot().mass(), centroidalInertia.block<3, 3>(0, 0));

  // (2) Update the virtual-contact model
  vcParams_.com = getRobot().com();
  // Use the value from the semi-axes-calculator 
  vcParams_.semiAxesVector << ssaPtr_->getSemiAxes()[0], ssaPtr_->getSemiAxes()[1], ssaPtr_->getSemiAxes()[2];
  // Use the impact body translation
  sva::PTransformd X_0_ee = getRobot().bodyPosW(getParams().iBodyName);
  vcParams_.eePosition = X_0_ee.translation(); 

  vcParams_.impactNormal = impactNormal;
  vcParams_.debug = false;
  vcPtr_->update(vcParams_);

  Eigen::Vector3d vc = vcPtr_->getVirtualContactPoint();

  // Update the twoDim model
  piParams_.contactPoint = ? 
  updatePlanarImpactParams_(
}

} // End of NameSpace
