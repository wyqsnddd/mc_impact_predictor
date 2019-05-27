#pragma once

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/Robots.h>

#include "mi_osd.h"

struct impulseValues
{

  bool inContact;
  Eigen::VectorXd deltaV;
  Eigen::VectorXd impulseForce;
  Eigen::VectorXd deltaTau;
  Eigen::VectorXd deltaQDot;

  void setContact()
  {
    inContact = true;
  }
  const bool contact()
  {
    return inContact;
  }

  void ini(const int & dim, const int & dof)
  {
    deltaV.resize(dim);
    impulseForce.resize(dim);
    deltaQDot.resize(dof);
    deltaTau.resize(dof);
    inContact = false;
  }

  void reset()
  {
    deltaV.setZero();
    impulseForce.setZero();
    deltaTau.setZero();
    deltaQDot.setZero();
  }
};
struct impactDataCache
{
  //
  Eigen::VectorXd qVelJump;
  Eigen::VectorXd tauJump;
  Eigen::MatrixXd jacobianDeltaAlpha;
  Eigen::MatrixXd jacobianDeltaTau;
  std::map<std::string, impulseValues> grfContainer;
  void reset()
  {
    qVelJump.setZero();
    tauJump.setZero();

    jacobianDeltaAlpha.setZero();
    jacobianDeltaTau.setZero();

    for(auto it = grfContainer.begin(); it != grfContainer.end(); ++it)
    {
      it->second.reset();
    }

  } // end of reset
  void ini(const int & dim, const int & dof)
  {
    qVelJump.resize(dof);
    tauJump.resize(dof);
    jacobianDeltaAlpha.resize(dof, dof);
    jacobianDeltaTau.resize(dof, dof);
    for(auto it = grfContainer.begin(); it != grfContainer.end(); ++it)
    {
      it->second.ini(dim, dof);
      it->second.reset();
    }

    reset();
  }
};

class mi_impactPredictor
{
public:
  mi_impactPredictor(mc_rbdyn::Robot & robot,
                     // std::shared_ptr<rbd::ForwardDynamics> & fdPtr,
                     std::string impactBodyName,
                     bool linearJacobian,
                     double impactDuration,
                     double coeRes = 0.8);

  ~mi_impactPredictor() {}
  void run(const Eigen::Vector3d & surfaceNormal);

  bool addEndeffector(std::string eeName);

  const Eigen::VectorXd & getTauJump() const
  {
    return cache_.tauJump;
  }
  const Eigen::VectorXd & getJointVelocityJump() const
  {
    return cache_.qVelJump;
  }
  /*
  const std::map<std::string, impulseValues>::iterator  nameToPointer(const std::string & eeName) const
  {

    const auto ee = cache_.grfContainer.find(eeName);
    if(ee != (cache_.grfContainer.end()))
    {
      return ee;
    }
    else
    {
      std::cout << "nameToPointer: " << eeName << std::endl;
      throw std::runtime_error("Predictor: required bodyname not found");
    }
  }
  */
  const Eigen::VectorXd & getEeVelocityJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaV;
  }

  const Eigen::VectorXd & getEeVelocityJump()
  {
    //  std::cout<<"calling eevel jump: "<<getImpactBody_()<<std::endl;
    const auto ee = cache_.grfContainer.find(getImpactBody_());

    // std::cout<<"Size of cache_.grfContainer is: "<<cache_.grfContainer.size()<<std::endl;
    /*
   if(ee==cache_.grfContainer.end()) {
    throw std::runtime_error("ee name is not found");
   }
    std::cout<<"impact body name is: "<<ee->first<<std::endl;
*/
    return ee->second.deltaV;
    // return Eigen::Vector3d::Zero();
  }
  const Eigen::VectorXd & getBranchJointVelJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaQDot;
  }
  const Eigen::VectorXd & getBranchTauJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaTau;
  }

  void resetDataStructure()
  {
    getOsd_()->resetDataStructure();
    cache_.reset();
  }
  void initializeDataStructure(int numEE)
  {

    cache_.ini(getOsd_()->getJacobianDim(), getOsd_()->getDof());
    getOsd_()->initializeDataStructure(numEE);
  }
  const Eigen::VectorXd & getImpulsiveForce()
  {
    const auto & ee = cache_.grfContainer.find(getImpactBody_());
    return ee->second.impulseForce;
  }

  const Eigen::VectorXd & getImpulsiveForce(const std::string & eeName)
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    if(ee->second.contact() || (ee->first == getImpactBody_()))
    {
      return ee->second.impulseForce;
    }
    else
    {
      throw std::runtime_error(std::string("Predictor: '-") + eeName + std::string("- ' is not in contact."));
    }
  }
  mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }
  /*
   mc_rbdyn::Robot & getRobot()
   {
     return robot_;
   }
 */
  void setImpactBody(std::string & impactBodyName)
  {
    impactBodyName_ = impactBodyName;
  }

  void setContact(const std::string contactBodyName)
  {
    const auto & ee = cache_.grfContainer.find(contactBodyName);
    ee->second.setContact();
    std::cout << "setContact: " << contactBodyName << ee->second.contact() << std::endl;
  }
  // We need to specify the endeffectors that are in contact, e.g. two feet
  // std::vector<dart::dynamics::BodyNode *> contactEndEffectors;

  //----------------Public API for control ----------------------------/
  const Eigen::MatrixXd & getJacobianDeltaAlpha()
  {
    return cache_.jacobianDeltaAlpha;
  }
  const Eigen::MatrixXd & getJacobianDeltaTau()
  {
    return cache_.jacobianDeltaTau;
  }
  const double & getImpactDuration_() const
  {
    return impactDuration_;
  }
protected:
  mc_rbdyn::Robot & robot_;
  std::string impactBodyName_;
  bool linearJacobian_;
  double impactDuration_;
  double coeRes_;

  std::shared_ptr<mi_osd> osdPtr_;

  impactDataCache cache_;

  bool useLinearJacobian_() const
  {
    return linearJacobian_;
  }

  // impact end-effector
  // dart::dynamics::BodyNodePtr impactBodyPtr_;

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
  double getCoeRes_() const
  {
    return coeRes_;
  }
  

  //  void tempTestEe_();
};
