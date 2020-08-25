#include "ImpactDynamicsModelInterface.h"
namespace mc_impact 
{

ImpactDynamicsModel::ImpactDynamicsModel(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params) : simRobot_(simRobot), params_(params)
{
  std::cout<<green<<"ImpactDynamicsModel is created."<<reset<<std::endl;

}

TwoDimModelBridge::TwoDimModelBridge(const mc_rbdyn::Robot & simRobot,
		const TwoDimModelBridgeParams & params) : ImpactDynamicsModel(simRobot, params.modelParams), params_(params)
{

  // Initialize the semiaxes calculator
  Eigen::Matrix3d inertiaMatrix = Eigen::Matrix3d::Identity();
  ssaPtr_.reset(new FIDynamics::SolveSemiAxes(getRobot().mass(), inertiaMatrix));
  
  // Initialize the virtual-contact-point calculator
  vcPtr_.reset(new FIDynamics::VirtualContactPoint());

  // Initialize the two-dim-model 
  // Energetic coefficient of restitution
  piParams_.e = 0.2;
  // Coefficient of friction 
  piParams_.miu = 0.7;
  rotation_.setZero();
  rotationFull_.setIdentity();

  twoDimModelPtr_.reset(new FIDynamics::TwoDimModel(piParams_));
  std::cout<<green<<"TwoDimModelBridge is created."<<reset<<std::endl;
}

/*
void TwoDimModelBridge::update( )
{

  // Update with the initial impact normal
  update(getParams().inertial_surfaceNormal);
}
*/
void TwoDimModelBridge::update(const Eigen::Vector3d & impactNormal, const Eigen::Vector3d & impactLinearVel)
{

  //std::cout<<"Updating TwoDimModelBridge"<<std::endl;
  // (0) Compute the whole-body inertia and average velocity
  Eigen::Matrix6d centroidalInertia; 
  centroidalInertia.setIdentity();

  //Eigen::Vector6d cm;
  Eigen::Vector6d av;
  rbd::computeCentroidalInertia(getRobot().mb(), getRobot().mbc(), getRobot().com(), centroidalInertia, av);
  
  //std::cout<<"inertia updated"<<std::endl;
  // Assert that the average com velocity is equal to the com velocity
  /*
  assert(mc_impact::areSame(av(3), getRobot().comVelocity()(0)));
  assert(mc_impact::areSame(av(4), getRobot().comVelocity()(1)));
  assert(mc_impact::areSame(av(5), getRobot().comVelocity()(2)));
  */

  rAverageAngularVel_ = av.segment<3>(0);
  rAverageLinearVel_ = av.segment<3>(3);
  rCentroidalInertia_ = centroidalInertia.block<3, 3>(0, 0);

  //std::cout<<"The average angular velocity is: "<<rAverageAngularVel_.transpose()<<std::endl;
  //std::cout<<"The average linear velocity is: "<<getRobot().comVelocity().transpose()<<std::endl;

  // (1) Update the ssa model
  
  // Inertia should be the upper corner? check! 
  //std::cout<<green<< "The centroidal inertia is: " << std::endl << centroidalInertia<< std::endl;
  ssaPtr_->update(getRobot().mass(), rCentroidalInertia_);

  //std::cout<<"ssa updated"<<std::endl;
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

  std::cout<<"The com is: "<<vcParams_.com.transpose()<<std::endl;
  // Use the value from the semi-axes-calculator 
  vcParams_.semiAxesVector << ssaPtr_->getSemiAxes()[0], ssaPtr_->getSemiAxes()[1], ssaPtr_->getSemiAxes()[2];
  std::cout<<"The semiAxesVector is: "<<vcParams_.semiAxesVector.transpose()<<std::endl;

  // Use the impact body translation
  sva::PTransformd X_0_ee = getRobot().bodyPosW(getParams().iBodyName);
  vcParams_.eePosition = X_0_ee.translation(); 

  std::cout<<"The ee position is: "<<vcParams_.eePosition.transpose()<<std::endl;
  vcParams_.impactNormal = impactNormal;
  vcParams_.debug = false;

  //std::cout<<"vc parameters updated"<<std::endl;
  
  Eigen::Vector3d vc; 
  if(getTwoDimModelBridgeParams().useVirtualContact)
  {
    vcPtr_->update(vcParams_);
    vc = vcPtr_->getVirtualContactPoint();
  }
  else
  {
    vc = vcParams_.eePosition;
  }

  std::cout<<"The virtual point is: "<<vc.transpose()<<std::endl;

  std::cout<<"The impact normal is: "<<vcParams_.impactNormal.transpose()<<std::endl;
  // (3) Update the twoDim model 
  updatePiParams_(vcParams_.impactNormal, vc, impactLinearVel);

  twoDimModelPtr_->updateParams(getPlanarImpactParams());
  //std::cout<<"twoDimModelPtr_->updated params "<<std::endl;
  twoDimModelPtr_->update();

  //std::cout<<"twoDimModelPtr_->updated "<<std::endl;
  // (4) Convert the twoDim model solution back to 3D:
  planarSolutionTo3D_();

  //std::cout<<"converted solution to 3d"<<std::endl;
}

void TwoDimModelBridge::updatePiParams_(const Eigen::Vector3d & in, const Eigen::Vector3d vc, const Eigen::Vector3d & impactLinearVel)
{
  // (1) Update the normal and tangential unit vectors
   // Compute the angle
   rotationAngle_= atan2(in.z(), in.y());
   std::cout<<green<<"The rotation angle is: "<<rotationAngle_<<std::endl;
   // Update the 2*3 rotation matrix:
   rotation_(0,0) = 1.0;
   rotation_(1,1) = cos(rotationAngle_);
   rotationFull_(1,1) = cos(rotationAngle_);
   rotation_(1,2) = -sin(rotationAngle_);
   rotationFull_(1,2) = -sin(rotationAngle_);

   rotationFull_(2,1) = sin(rotationAngle_);
   rotationFull_(2,2) = cos(rotationAngle_);

   piParams_.nu = rotation_*in;
   //Eigen::Vector2d rotatedZ = rotation_*Eigen::Vector3d::UnitZ();
   piParams_.tu(0) = -piParams_.nu(1);
   piParams_.tu(1) = piParams_.nu(0);

   std::cout<<green<<"The nu is: "<<piParams_.nu.transpose()<<std::endl;
   std::cout<<green<<"The tu is: "<<piParams_.tu.transpose()<<std::endl;

  // (2) Contact Point: 
   piParams_.contactPoint = rotation_*vc;

   std::cout<<green<<"The rotated contact point is: "<<piParams_.contactPoint.transpose()<<std::endl;
  // (3) Parmams of the bat and the object:
   switch(getCase_())
   {
     case TwoDimModelCase::PushWall:
	     // Update the parameters using the Push-Wall assumptions
	     paramUpdatePushWall_(impactLinearVel);
	     break;
     default:
	  throw std::runtime_error("The assumptions are not set for the TwoDimModelBridge.");
   }
  
}

void TwoDimModelBridge::paramUpdatePushWall_(const Eigen::Vector3d & impactLinearVel)
{
  // Bat is supposed to be the robot:
  // (1) Robot
  piParams_.batParams.com = rotation_*getRobot().com(); 
  std::cout<<green<<"The rotated com is: "<<piParams_.batParams.com.transpose()<<std::endl;
  piParams_.batParams.mass = getRobot().mass(); 
  // Get the z-axis diagonal element: 
  piParams_.batParams.inertia = (rotationFull_*rCentroidalInertia_)(2,2); 
  std::cout<<green<<"The bat inertia is: "<<piParams_.batParams.inertia<<std::endl;
  // Rotate the com velocity 
  if(getTwoDimModelBridgeParams().useComVel)
  {
    piParams_.batParams.preVel = rotation_*rAverageLinearVel_;
    //piParams_.batParams.preVel = rotation_*getRobot().bodySensor("FloatingBase").linearVelocity();
  }
  else
  {
    piParams_.batParams.preVel = rotation_*impactLinearVel;
  }

  std::cout<<green<<"The bat preimpact vel is: "<<piParams_.batParams.preVel<<std::endl;
  // Get the z-axis average angular velocity:
  piParams_.batParams.preW = rAverageAngularVel_(2);
  std::cout<<green<<"The bat preimpact angular vel is: "<<piParams_.batParams.preW<<std::endl;

  piParams_.batParams.name = "robot";
  
  // Object is suppose to be the wall:
  
  // (2) Wall
  piParams_.objectParams.preVel << 0.0, 0.0;
  piParams_.objectParams.preW = 0.0;
  piParams_.objectParams.name = "wall";

  // mass and inertia of the wall are set to be infinite.
  //piParams_.objectParams.mass = std::numeric_limits<double>::infinity();
  piParams_.objectParams.mass = std::numeric_limits<double>::max();
  //piParams_.objectParams.inertia = std::numeric_limits<double>::infinity();
  piParams_.objectParams.inertia = std::numeric_limits<double>::max();

  
  // com of the wall is set to be the contact point such that r = cp - com == 0.
  // Suppose that piParams_.contactPoint is already set.
  piParams_.objectParams.com = piParams_.contactPoint;
}

void TwoDimModelBridge::planarSolutionTo3DPushWall_()
{
  // Convert the post-impact impulse:
  // The robot applies the impulse "I", thus it receives impulse "-I". 
  robotPostImpactStates_.impulse = - rotation_.transpose()*twoDimModelPtr_->getSolution().I_r;
  std::cout<<"I_nr is: "<< twoDimModelPtr_->getSolution().I_nr<<std::endl;
  std::cout<<"I_nc is: "<< twoDimModelPtr_->getSolution().I_nc<<std::endl;
  std::cout<<"I_r is: "<< twoDimModelPtr_->getSolution().I_r.transpose()<<std::endl;
  std::cout<<"The impulse is: "<< robotPostImpactStates_.impulse<<std::endl;
  // Convert the post-impact velocities:
  // robot:
  robotPostImpactStates_.linearVel = rotation_.transpose()*twoDimModelPtr_->getImpactBodies().first.postVel;
  robotPostImpactStates_.linearVelJump = robotPostImpactStates_.impulse/getRobot().mass();

  std::cout<<"The post-imapct linear velocity is: "<< twoDimModelPtr_->getImpactBodies().first.postVel<<std::endl;

  //  Compute: wJump = (1 / getParams().inertia) * cross2(r, I_r);
  Eigen::Vector3d rb = vcPtr_->getVirtualContactPoint() - vcParams_.com;
  robotPostImpactStates_.anguleVelJump = rCentroidalInertia_.inverse()*rb.cross(robotPostImpactStates_.impulse);

  robotPostImpactStates_.anguleVel = rAverageAngularVel_ + robotPostImpactStates_.anguleVelJump;

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


void TwoDimModelBridge::logImpulseEstimations()
{
  const std::string & bridgeName = getTwoDimModelBridgeParams().name; 

  logEntries_.emplace_back(bridgeName + "_"+ "rotationAngle");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return rotationAngle_; });

  logEntries_.emplace_back(bridgeName + "_"+ "robot_postImpact_linearVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().linearVel; });

  logEntries_.emplace_back(bridgeName + "_"+ "robot_linearVelJump");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().linearVelJump; });


  logEntries_.emplace_back(bridgeName + "_"+ "robot_postImpact_angularVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().anguleVel; });

  logEntries_.emplace_back(bridgeName + "_"+ "robot_postImpact_angularVelJump");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().anguleVelJump; });


  logEntries_.emplace_back(bridgeName + "_"+ "robot_postImpact_impulse");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().impulse; });

  std::cout<<green<<"The impulse is: "<<getRobotPostImpactStates().impulse<<std::endl;

  logEntries_.emplace_back(bridgeName + "_"+ "average_angularVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return rAverageAngularVel_; });

  logEntries_.emplace_back(bridgeName + "_"+ "average_linearVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return rAverageLinearVel_; });


  logEntries_.emplace_back(bridgeName + "_"+ "virtualContactPoint");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return vcPtr_->getVirtualContactPoint(); });

  logEntries_.emplace_back(bridgeName + "_"+ "useVirtualContactPoint");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { 
		 if(getTwoDimModelBridgeParams().useVirtualContact) 
		  return true;
		 else 
		  return false; 
		  });

  logEntries_.emplace_back(bridgeName + "_"+ "useComVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { 
		 if(getTwoDimModelBridgeParams().useComVel) 
		  return true;
		 else 
		  return false; 
		  });

  logEntries_.emplace_back(bridgeName + "_"+ "impactEventSequence");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { 
		  switch(twoDimModelPtr_->getSolution().sequence){
		     case FIDynamics::EventSequence::not_specified:
		         return 0; 
		     case FIDynamics::EventSequence::C_R_Process:
			 return 1;
		     case FIDynamics::EventSequence::C_S_R_Process:
			 return 2;
		     case FIDynamics::EventSequence::S_C_R_Process:
			 return 3;
		     case FIDynamics::EventSequence::SpecialSolution:
			 return 4;
		  default:  
		      throw std::runtime_error("impact event sequence is not defined");
		  }
		 });
  /*
  logEntries_.emplace_back(bridgeName + "_"+ "object_postImpact_linearVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getObjectPostImpactStates().linearVel; });

  logEntries_.emplace_back(bridgeName + "_"+ "object_postImpact_angularVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getObjectPostImpactStates().anguleVel; });

  logEntries_.emplace_back(bridgeName + "_"+ "object_postImpact_impulse");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getObjectPostImpactStates().impulse; });
*/


}
void TwoDimModelBridge::removeImpulseEstimations()
{
  assert(getHostCtl_() not nullptr);

  for(auto & name : logEntries_)
  {
    getHostCtl_()->logger().removeLogEntry(name); 
  }
}
} // End of NameSpace
