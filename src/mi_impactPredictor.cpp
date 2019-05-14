#include "mi_impactPredictor.h"

mi_impactPredictor::mi_impactPredictor(const dart::dynamics::SkeletonPtr & robotPtr,
				       // const mc_rbdyn::Robot & robot,
                                       const std::string & impactBodyName,
                                       bool linearJacobian,
                                       // const dart::dynamics::BodyNodePtr impactBodyPtr,
                                       double impactDuration,
                                       double coeRes)
: robotPtr_(robotPtr), linearJacobian_(linearJacobian), impactDuration_(impactDuration), coeRes_(coeRes)
{

  std::cout << "The impact predictor constuctor is started." << std::endl;
  setImpactBody(impactBodyName);
  //osdPtr_ = std::make_shared<mi_osd>(robotPtr, getRobot(), useLinearJacobian_());
  osdPtr_ = std::make_shared<mi_osd>( getDartRobot(), useLinearJacobian_());
  std::cout << "The impact predictor constuctor is finished." << std::endl;
  std::cout << "The impact duration is: " << getImpactDuration_() << ", the coeres is: " << getCoeRes_() << std::endl;
  cache_.eeVelJump.setZero();
  cache_.qVelJump.setZero();
  cache_.eeImpulse.setZero();
  // Initialize the GRF container: <delta-Vel, delta-Impulse>
  // std::string l_foot_name("l_ankle");
  // std::string r_foot_name("r_ankle");
  if(useLinearJacobian_())
  {
    cache_.grfContainer["l_ankle"] =
        std::make_pair<Eigen::VectorXd, Eigen::VectorXd>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    cache_.grfContainer["r_ankle"] =
        std::make_pair<Eigen::VectorXd, Eigen::VectorXd>(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
  }
  else
  {
    cache_.grfContainer["l_ankle"] =
        std::make_pair<Eigen::Vector6d, Eigen::Vector6d>(Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero());
    cache_.grfContainer["r_ankle"] =
        std::make_pair<Eigen::Vector6d, Eigen::Vector6d>(Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero());
  }
}

void mi_impactPredictor::run()
{
  // Update the equations of motions
  std::cout << "mi_impactPredictor::update() is called. " << std::endl;
  osdPtr_->update();
  // Update the data in the cache
  // * Update the end-effector velocityJump
  std::cout << "OSD updated. " << std::endl;
  //Eigen::VectorXd alpha = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  //Eigen::VectorXd alphaD = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alphaD);
  std::cout << "q_d is: " << getDartRobot()->getVelocities()<< std::endl;
  std::cout << "The velocity of " << getImpactBody_() << " is: " << getDartRobot()->getBodyNode(getImpactBody_())->getSpatialVelocity() << std::endl;
  // std::cout<<"q_dd is: "<<alphaD.transpose()<<std::endl;
  // std::cout<<" The Jacobian is: "<<getOsd_()->getJacobian(getImpactBody_())<<std::endl;

  cache_.eeVelJump =
      -(getCoeRes_() + 1) * getOsd_()->getJacobian(getImpactBody_()) * 
      (getDartRobot()->getVelocities() + getDartRobot()->getAccelerations()* getImpactDuration_());

  // std::cout<<"The impact body is: "<<getImpactBody_()<<std::endl;
  // * Update the joint velocity jump
  cache_.qVelJump = getOsd_()->getDcJacobianInv(getImpactBody_()) * cache_.eeVelJump;
  std::cout << "The predicted impact body velocity jump is: " << std::endl << cache_.eeVelJump.transpose() << std::endl;
  // * Update the impulsive force
  /*
  cache_.eeImpulse = (1 / getImpactDuration_())
    *getOsd_()->getEffectiveLambdaMatrix(getImpactBody_())
    * getOsd_()->getInvMassMatrix().transpose()
    * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_()).transpose();
   */
  cache_.eeImpulse = (1 / getImpactDuration_()) 
	  * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_())
	  * getOsd_()->getDcJacobianInv(getImpactBody_()) 
	  * cache_.eeVelJump;

  // * Update the impulsive force of end-effectors with established contact

  // std::cout<<"The predicted joint velocity jump is: "<<std::endl<<cache_.qVelJump.transpose()<<std::endl;
  std::cout << "The predicted end-effector impulsive force is: " << std::endl
            << cache_.eeImpulse.transpose() << std::endl;

  // Update the ground reaction forces:
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {

    // End-effector velocity jump:
    it->second.first = getOsd_()->getJacobian(it->first) 
	    * cache_.qVelJump;

    // End-effector reaction force:
    it->second.second = (1 / getImpactDuration_()) 
	    * getOsd_()->getEffectiveLambdaMatrix(it->first)
	    * getOsd_()->getDcJacobianInv(getImpactBody_()) 
	    * cache_.eeVelJump;

    std::cout << "The predicted GRF impulsive force of " << it->first << " is: " << std::endl
              << it->second.second.transpose() << ", velocity jump is: " << it->second.first.transpose() << std::endl;
  }
  tempTest_();
}



void mi_impactPredictor::tempTest_(){

// Test the end-effector induced ground reaction forces. 

// Note that we need to deduct the gravity force. 
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
/*	
    // Read the acceleration of the impact body: 
    auto tempImpactBodyAcceleration = 
	    getDartRobot()->getBodyNode(it->first)
	    getRobot().mbc().bodyPosW[getRobot().mb().bodyIndexByName(it->first)]
	    *getRobot().mbc().bodyAccB
	    [
	    getRobot().mb().bodyIndexByName(it->first)
	    ];
    // End-effector reaction force:
    auto tempGRF = (1 / getImpactDuration_()) 
	    *getOsd_()->getEffectiveLambdaMatrix(it->first)
	    * getOsd_()->getDcJacobianInv(getImpactBody_())
	    * tempImpactBodyAcceleration.linear();
    std::cout<<"The impact body acceleration is: "<<tempImpactBodyAcceleration.linear()<<std::endl;
    std::cout << "The predicted GRF force due to impact-body movement of: "
	    << it->first << " is: " << std::endl
              << tempGRF << std::endl;
	      */
  }
 

}


