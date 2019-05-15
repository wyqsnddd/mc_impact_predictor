#pragma once

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/Robots.h>

//#include "mc_dart_controller.h"
#include "mi_osd.h"
// #include <dart/constraint/constraint.hpp>
// #include <dart/dynamics/dynamics.hpp>


struct impulseValues{
	Eigen::VectorXd deltaV;
	Eigen::VectorXd impulseForce;
	Eigen::VectorXd accForce;
};
struct impactDataCache
{
  //
  Eigen::VectorXd eeVelJump;
  Eigen::VectorXd qVelJump;
  Eigen::VectorXd tauJump;
  Eigen::VectorXd eeImpulse;
  /// <end-effector Name, <delta-V, delta-F, F-due-to-ee-acc>>
  std::map<std::string, impulseValues > grfContainer;
  //std::map<std::string, std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> > grfContainer;
};

class mi_impactPredictor
{
public:
  mi_impactPredictor(// const dart::dynamics::SkeletonPtr & robotPtr,
		     mc_rbdyn::Robot & robot,
                     const std::string & impactBodyName,
                     bool linearJacobian,
                     double impactDuration,
                     double coeRes = 0.8);

  ~mi_impactPredictor() {}

  void run();

  Eigen::VectorXd getEeVelocityJump() const
  {
    return cache_.eeVelJump;
  }

  Eigen::VectorXd getEeVelocityJump(const std::string & eeName) const
  {
    auto ee = cache_.grfContainer.find(eeName);
    //return ee->second.first;
    return ee->second.deltaV;
  }

  Eigen::VectorXd getJointVelocityJump() const
  {
    return cache_.qVelJump;
  }
  Eigen::VectorXd getImpulsiveForce() const
  {
    return cache_.eeImpulse;
  }

  Eigen::VectorXd getImpulsiveForce(const std::string & eeName) const
  {
    auto ee = cache_.grfContainer.find(eeName);
    //return ee->second.second;
    return ee->second.impulseForce;
  }

  Eigen::VectorXd getEeAccForce(const std::string & eeName) const
  {
    auto ee = cache_.grfContainer.find(eeName);
    return ee->second.accForce;
  }

  const mc_rbdyn::Robot & getRobot() const
  {
    return robot_;
  }
  mc_rbdyn::Robot & getRobot() 
  {
    return robot_;
  }
  void setImpactBody(std::string impactBodyName)
  {
    impactBodyName_ = impactBodyName;
  }
  // We need to specify the endeffectors that are in contact, e.g. two feet
  // std::vector<dart::dynamics::BodyNode *> contactEndEffectors;

protected:
  std::shared_ptr<mi_osd> osdPtr_;

  mc_rbdyn::Robot & robot_;
  bool useLinearJacobian_() const
  {
    return linearJacobian_;
  }
  bool linearJacobian_;

  // impact end-effector
  // dart::dynamics::BodyNodePtr impactBodyPtr_;
  std::string impactBodyName_;

  // The returned pointer is not supposed to be changed.
  const std::string getImpactBody_()
  {
    return impactBodyName_;
  }

  // The returned pointer is not supposed to be changed.
  const std::shared_ptr<mi_osd> getOsd_() const
  {
    return osdPtr_;
  }

  // mc_rbdyn::Robot & loadRobot_(mc_rbdyn::RobotModulePtr rm, const std::string & name);
  double impactDuration_;
  double coeRes_;
  double getCoeRes_() const
  {
    return coeRes_;
  }
  double getImpactDuration_() const
  {
    return impactDuration_;
  }

  impactDataCache cache_;
  void tempTest_();
};
