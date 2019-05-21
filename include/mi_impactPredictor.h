#pragma once

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/Robots.h>

#include "mi_osd.h"

struct impulseValues
{
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
  Eigen::VectorXd newLeeImpulse;
  Eigen::VectorXd newReeImpulse;
  Eigen::VectorXd new_eeLeeImpulse;
  Eigen::VectorXd new_eeReeImpulse;
  /// <end-effector Name, <delta-V, delta-F, F-due-to-ee-acc>>
  std::map<std::string, impulseValues> grfContainer;
  // std::map<std::string, std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> > grfContainer;
};

class mi_impactPredictor
{
public:
  mi_impactPredictor(mc_rbdyn::Robot & robot,
                     const std::string & impactBodyName,
                     bool linearJacobian,
                     double impactDuration,
                     double coeRes = 0.8);

  ~mi_impactPredictor() {}
  //void initializeDataStructure();
  void resetDataStructure();
  void run();
  bool addEndeffector(const std::string & eeName)
  {

    unsigned eeNum = static_cast<unsigned>(cache_.grfContainer.size());

    if(useLinearJacobian_())
    {
      cache_.grfContainer[eeName] = {Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()};
    }
    else
    {
      cache_.grfContainer[eeName] = {Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero()};
    }

    if(!getOsd_()->addEndeffector(eeName))
    {
      throw std::runtime_error("OSd failed to add endeffector!");
    }

    unsigned eeNumNew = static_cast<unsigned>(cache_.grfContainer.size());
    if((eeNum == eeNumNew - 1) && (eeNumNew == static_cast<unsigned>(getOsd_()->getEeNum())))
    {
      return true;
    }
    else
    {
      std::cout << "The ee container start with " << eeNum << ", now it has " << eeNumNew
                << " endeffectors. The osd has " << getOsd_()->getEeNum() << " end-effectors. " << std::endl;
      return false;
    }
  }
  const Eigen::VectorXd & getEeVelocityJump() const
  {
    return cache_.eeVelJump;
  }
  const Eigen::VectorXd & getTauJump() const
  {
    return cache_.tauJump;
  }
  const Eigen::VectorXd & getEeVelocityJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaV;
  }

  const Eigen::VectorXd & getJointVelocityJump() const
  {
    return cache_.qVelJump;
  }
  const Eigen::VectorXd & getImpulsiveForce() const
  {
    return cache_.eeImpulse;
  }
  const Eigen::VectorXd & getNewLeeImpulsiveForce() const
  {
    return cache_.newLeeImpulse;
  }
  const Eigen::VectorXd & getNewReeImpulsiveForce() const
  {
    return cache_.newReeImpulse;
  }
  const Eigen::VectorXd & getNewEeLeeImpulsiveForce() const
  {
    return cache_.new_eeLeeImpulse;
  }
  const Eigen::VectorXd & getNewEeReeImpulsiveForce() const
  {
    return cache_.new_eeReeImpulse;
  }
  const Eigen::VectorXd & getImpulsiveForce(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.impulseForce;
  }

  const Eigen::VectorXd & getEeAccForce(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
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
  void tempTestBody_();
  void tempTestEe_();
  void tempTestAcc_();
  void tempTestAccEe_();
};
