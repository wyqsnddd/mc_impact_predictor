#pragma once

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/Robots.h>

#include "mi_osd.h"

struct impulseValues
{

  bool inContact;
  // Eigen::Vector3d deltaCoP;
  Eigen::VectorXd deltaV;
  Eigen::VectorXd impulseForce;
  /// This is the equivalent impulsive wrench at the COM
  sva::ForceVecd impulseForceCOM;
  sva::ForceVecd impulseForceCOP;
  // Eigen::VectorXd accForce;
  Eigen::VectorXd deltaTau;
  Eigen::VectorXd deltaQDot;
  Eigen::MatrixXd jacobianDeltaF;

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
    // accForce.resize(dim);
    deltaQDot.resize(dof);
    deltaTau.resize(dof);
    inContact = false;
  }

  void reset()
  {
    deltaV.setZero();
    impulseForce.setZero();
    impulseForceCOM.vector().setZero();
    // deltaCoP.setZero();
    // accForce.setZero();
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
  std::vector<std::string> contactEndeffectors;
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
                     double coeFrictionDeduction,
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
  const Eigen::VectorXd & getJointVelocityJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaQDot;
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
  /*
  const Eigen::Vector3d & getDeltaCoP(const std::string & eeName)
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaCoP;
  }
*/
  /*
  Eigen::VectorXd getOsdForce(const std::string & eeName)
  {
    return getOsd_()->getOsdForce(eeName);
  }
  */
  /*
   const Eigen::VectorXd & getAccForce(const std::string & eeName)
   {
     const auto & ee = cache_.grfContainer.find(eeName);
     return ee->second.accForce;
   }
   */
  /*
  Eigen::VectorXd getQPForce(const std::string & eeName)
  {
    std::cout<<"The dc jacobian size is: "<< getOsd_()->getDcJacobianInv(eeName).transpose().size()<<std::endl;
    std::cout<<"The joint torque size is: "<< rbd::dofToVector(getRobot().mb(),
  getRobot().mbc().jointTorque).size()<<std::endl; return getOsd_()->getDcJacobianInv(eeName).transpose()
      *rbd::dofToVector(getRobot().mb(), getRobot().mbc().jointTorque);
  }
  */
  const Eigen::VectorXd & getImpulsiveForce()
  {
    const auto & ee = cache_.grfContainer.find(getImpactBody_());
    return ee->second.impulseForce;
  }
  const sva::ForceVecd & getImpulsiveForceCOM()
  {
    /*
    sva::PTransformd X_ee_CoM = sva::PTransformd(getRobot().com())*getRobot().bodyPosW(getImpactBody_()).inv();

    return X_ee_CoM.dualMul(sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce()));
    */
    const auto & ee = cache_.grfContainer.find(getImpactBody_());
    return ee->second.impulseForceCOM;
  }
  const sva::ForceVecd & getImpulsiveForceCOP(const std::string & eeName)
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    if(ee->second.contact())
    {
      return ee->second.impulseForceCOP;
    }
    else
    {
      throw std::runtime_error(std::string("Predictor: '-") + eeName + std::string("- ' is not in contact."));
    }
  }
  const sva::ForceVecd & getImpulsiveForceCOM(const std::string & eeName)
  {
    /*
    sva::PTransformd X_ee_CoM = sva::PTransformd(getRobot().com())*getRobot().bodyPosW(eeName).inv();
    return X_ee_CoM.dualMul(sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce(eeName)));
    */
    const auto & ee = cache_.grfContainer.find(eeName);
    if(ee->second.contact() || (ee->first == getImpactBody_()))
    {
      return ee->second.impulseForceCOM;
    }
    else
    {
      throw std::runtime_error(std::string("Predictor: '-") + eeName + std::string("- ' is not in contact."));
    }
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

  void setContact(std::string contactBodyName)
  {
    const auto & ee = cache_.grfContainer.find(contactBodyName);
    if(ee != (cache_.grfContainer.end()))
    {
      ee->second.setContact();
      cache_.contactEndeffectors.push_back(contactBodyName);
      std::cout << "setContact: " << contactBodyName << ee->second.contact() << std::endl;
    }
    else
    {
      std::cout << "setContact: " << contactBodyName << std::endl;
      throw std::runtime_error(std::string("setContact: '-") + contactBodyName
                               + std::string("- ' is not in the container."));
    }
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
  const Eigen::MatrixXd & getJacobianDeltaF(const std::string & eeName)
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.jacobianDeltaF;
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

  double coeFrictionDeduction_;
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
  const double getCoeRes_() const
  {
    return coeRes_;
  }
  const double getCoeFricDe_() const
  {
    return coeFrictionDeduction_;
  }
  //  void tempTestEe_();
};
