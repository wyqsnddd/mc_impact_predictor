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
  rotation_.setZero();
  rotationFull_.setIdentity();

  twoDimModelPtr_.reset(new FIDynamics::TwoDimModel(piParams_));
  std::cout<<green<<"TwoDimModelBridge is created."<<reset<<std::endl;
}

void TwoDimModelBridge::update( )
{

  // Update with the initial impact normal
  update(getParams().inertial_surfaceNormal);
}
void TwoDimModelBridge::update(const Eigen::Vector3d & impactNormal)
{

  // (0) Compute the whole-body inertia and average velocity
  Eigen::Matrix6d centroidalInertia; 
  centroidalInertia.setIdentity();

  Eigen::Vector6d cm;
  Eigen::Vector6d av;
  rbd::computeCentroidalInertiaAndVelocity(getRobot().mb(), getRobot().mbc(), getRobot().com(), centroidalInertia, cm, av);
  
  // Assert that the average com velocity is equal to the com velocity
  assert(mc_impact::areSame(av(3), getRobot().comVelocity()(0)));
  assert(mc_impact::areSame(av(4), getRobot().comVelocity()(1)));
  assert(mc_impact::areSame(av(5), getRobot().comVelocity()(2)));

  rAverageAngularVel_ = av.segment<3>(0);
  rCentroidalInertia_ = centroidalInertia.block<3, 3>(0, 0);

  // (1) Update the ssa model
  
  // Inertia should be the upper corner? check! 
  //std::cout<<green<< "The centroidal inertia is: " << std::endl << centroidalInertia<< std::endl;
  ssaPtr_->update(getRobot().mass(), rCentroidalInertia_);

// clang-format off
// This is a sample centroidal inertia:
/*
   6.50171    0.0572213     0.282871   1.1147e-17  7.68482e-16 -3.81639e-16
   0.0572213      6.06731       0.1157 -2.58474e-16  3.37024e-17 -3.95517e-16
   0.282871       0.1157      1.81269  1.94289e-16  1.78677e-16 -1.89735e-17
   1.8384e-17  2.13371e-16   3.1225e-17      39.0255  9.96488e-17 -1.82363e-16
   2.75821e-16  5.28301e-17  1.17961e-16  4.28319e-17      39.0255  1.39119e-16
  -5.89806e-16 -4.02456e-16 -1.10453e-18 -6.34367e-16  2.32848e-16      39.0255
*/
// clang-format on
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

  // (3) Update the twoDim model 
  updatePiParams_(vcParams_.impactNormal, vc);
  twoDimModelPtr_->update();

  // (4) Convert the twoDim model solution back to 3D:
  planarSolutionTo3D_();

}

void TwoDimModelBridge::updatePiParams_(const Eigen::Vector3d & in, const Eigen::Vector3d vc)
{
  // (1) Update the normal and tangential unit vectors
   // Compute the angle
   double angle = atan2(in.z(), in.y());
   // Update the 2*3 rotation matrix:
   rotation_(0,0) = 1.0;
   rotation_(1,1) = cos(angle);
   rotationFull_(1,1) = cos(angle);
   rotation_(1,2) = -sin(angle);
   rotationFull_(1,2) = -sin(angle);

   rotationFull_(2,1) = sin(angle);
   rotationFull_(2,2) = cos(angle);

   piParams_.nu = rotation_*in;
   //Eigen::Vector2d rotatedZ = rotation_*Eigen::Vector3d::UnitZ();
   piParams_.tu(0) = -piParams_.nu(1);
   piParams_.tu(1) = piParams_.nu(0);

  // (2) Contact Point: 
   piParams_.contactPoint = rotation_*vc;

  // (3) Parmams of the bat and the object:
   switch(getCase_())
   {
     case TwoDimModelCase::PushWall:
	     // Update the parameters using the Push-Wall assumptions
	     paramUpdatePushWall_();
	     break;
     default:
	  throw std::runtime_error("The assumptions are not set for the TwoDimModelBridge.");
   }
  
}

void TwoDimModelBridge::paramUpdatePushWall_()
{
  // Bat is supposed to be the robot:
  // (1) Robot
  piParams_.batParams.com = rotation_*getRobot().com(); 
  piParams_.batParams.mass = getRobot().mass(); 
  // Get the z-axis diagonal element: 
  piParams_.batParams.inertia = (rotationFull_*rCentroidalInertia_)(2,2); 
  // Rotate the com velocity 
  piParams_.batParams.preVel = rotation_*getRobot().comVelocity();
  // Get the z-axis average angular velocity:
  piParams_.batParams.preW = rAverageAngularVel_(2);

  piParams_.batParams.name = "robot";
  
  // Object is suppose to be the wall:
  
  // (2) Wall
  piParams_.objectParams.preVel << 0.0, 0.0;
  piParams_.objectParams.preW = 0.0;
  piParams_.objectParams.name = "wall";

  // mass and inertia of the wall are set to be infinite.
  piParams_.objectParams.mass = std::numeric_limits<double>::infinity();
  piParams_.objectParams.inertia = std::numeric_limits<double>::infinity();

  
  // com of the wall is set to be the contact point such that r = cp - com == 0.
  // Suppose that piParams_.contactPoint is already set.
  piParams_.objectParams.com = piParams_.contactPoint;
}

void TwoDimModelBridge::planarSolutionTo3DPushWall_()
{
  // Convert the post-impact impulse:
  // The robot applies the impulse "I", thus it receives impulse "-I". 
  robotPostImpactStates_.impulse = - rotation_.transpose()*twoDimModelPtr_->getSolution().I_r;

  // Convert the post-impact velocities:
  // robot:
  robotPostImpactStates_.linearVel = rotation_.transpose()*twoDimModelPtr_->getImpactBodies().first.postVel;
  robotPostImpactStates_.anguleVel = rAverageAngularVel_;
  robotPostImpactStates_.anguleVel(2) += twoDimModelPtr_->getImpactBodies().first.postW;

}
void TwoDimModelBridge::planarSolutionTo3D_()
{
  switch(getCase_())
   {
     case TwoDimModelCase::PushWall:
	     // Convert the solution using the Push-Wall assumptions
	     planarSolutionTo3DPushWall_();
	     break;
     default:
	  throw std::runtime_error("The assumptions are not set for the TwoDimModelBridge.");
   }
}

const PostImpactStates & TwoDimModelBridge::getObjectPostImpactStates()
{
   switch(getCase_())
   {
     case TwoDimModelCase::PushWall:
	     // In this case the object(wall)  is supposed to be stationary.
	     throw std::runtime_error("In the PushWall case, the wall is stationary. Thus there is no need to check its post-impact states.");
     default:
	  return objectPostImpactStates_; 
   }
 }


const PostImpactStates & ImpactDynamicsModel::getRobotPostImpactStates()
{
   return robotPostImpactStates_; 
}


const PostImpactStates & ImpactDynamicsModel::getObjectPostImpactStates()
{
   return objectPostImpactStates_; 
}

} // End of NameSpace
