#pragma once

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/Robots.h>

#include "mi_osd.h"

struct impulseValues
/** \brief Impulsive values associated with one end-effector
 */
{
  bool inContact; ///< If this end-effector is in contact

  Eigen::VectorXd deltaV; ///< End-effector velocity jump \f$ \delta \v_k \f$

  Eigen::VectorXd impulseForce; ///< End-effector velocity jump \f$ \bar{f}_k \f$

  sva::ForceVecd impulseForceCOM; ///< Equivalent force of \f$ \bar{f}_k \f$ at the CoM

  sva::ForceVecd impulseForceCOP;
  // Eigen::VectorXd accForce;
  Eigen::VectorXd deltaTau; ///< Impact induced joint torque jump \f$ \delta \tau_k \f$ of the \f$k\f$th kinematic
                            ///< branch associated with the end-effector */

  Eigen::VectorXd deltaQDot; ///< Impact induced joint velocity jump \f$ \delta \dot{q}_k \f$ of the \f$k\f$th kinematic
                             ///< branch associated with the end-effector

  Eigen::MatrixXd
      jacobianDeltaF; ///<   We use this Jacobian to calculate the impact-induced impulsive force: \f$ \bar{f}_m =
                      ///<   \frac{1}{\delta t}\mathcal{J}_{f}(\dot{q} + \ddot{q}\Delta t) \f$ where we have the
                      ///<   definition:   \f$  \mathcal{J}_{f} = \frac{1}{\delta t }  \sum^{m}_{i=1}\Lambda_{ki}
                      ///<   \Lambda^{-1}_{km} \tilde{\Lambda}_m P_{\delta}J_m. \f$

  /*!
   * Set this end-effector in contact
   */
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
/** \brief Impulsive values associated with the robot
 */
{
  Eigen::VectorXd qVelJump; ///< Impact induced joint velocity jump \f$ \delta \dot{q}\f$ of all the kinematic branches

  Eigen::VectorXd tauJump; ///< Impact induced joint torque jump \f$ \delta \tau \f$ of all the kinematic branches.

  Eigen::MatrixXd
      jacobianDeltaAlpha; ///<     We use this Jacobian to calculate the impact-induced joint velocity jump: \f$ \delta
                          ///<     \dot{q} = \frac{1}{\delta t}\mathcal{J}_{\delta \dot{q}}(\dot{q} + \ddot{q}\Delta t)
                          ///<     \f$ where we have the definition: \f$ \mathcal{J}_{\delta \dot{q}} = ( \bar{J}_m  +
                          ///<     \frac{1}{\delta t} (\sum_{i \neq m}  \bar{J}_i \Lambda^{-1}_{im}
                          ///<     )\tilde{\Lambda}_{m}) P_{\delta}J_m \f$

  Eigen::MatrixXd
      jacobianDeltaTau; ///< We use this Jacobian to calculate the impact-induced joint torque jump: \f$ \delta \tau =
                        ///< \frac{1}{\delta t}\mathcal{J}_{\delta \tau }(\dot{q} + \ddot{q}\Delta t) \f$ where we have
                        ///< the definition: \f$  \mathcal{J}_{\delta \tau} = (J^\top_m  +  \frac{1}{\delta t }
                        ///< \sum_{k\in \mathcal{I}_c} J^\top_k\sum^{m}_{i=1}\Lambda_{ki}  \Lambda^{-1}_{km}
                        ///< )\tilde{\Lambda}_m P_{\delta}J_m \f$

  std::map<std::string, impulseValues>
      grfContainer; ///< a container that helps to find the impulsive values associated with an end-effector.
  std::vector<std::string> contactEndeffectors; ///< end-effectors with established contact.
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
/** \brief impact-induced state jumps calculation when one end-effector is expecting an impact.
 */
{
public:
  /*!
      \param robot reference to the robot model used by the QP
      \param osdPtr reference to the shared pointer of the Operational space dynamics
      \param impactBodyName name of the end-effector that expects an impact
      \param coeFrictionDeduction tangential velocity deduction coeffecient at the impact.
      \param coeRes coefficient of restitution
      */
  mi_impactPredictor(mc_rbdyn::Robot & robot,
                     const std::shared_ptr<mi_osd> & osdPtr,
                     std::string impactBodyName,
                     bool linearJacobian,
                     double impactDuration,
                     double coeFrictionDeduction,
                     double coeRes = 0.8);

  ~mi_impactPredictor() {}
  /*!
  \param surfaceNormal normal direction of the expected contact surface.
  */
  void run(const Eigen::Vector3d & surfaceNormal);

  /*!
    \param add the end-effector to the operational space model
    \return operation status
    */
  bool addEndeffector(std::string eeName);

  inline const Eigen::VectorXd & getTauJump() const
  {
    return cache_.tauJump;
  }
  inline const Eigen::VectorXd & getJointVelocityJump() const
  {
    return cache_.qVelJump;
  }
  inline const Eigen::VectorXd & getJointVelocityJump(const std::string & eeName) const
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
  inline const Eigen::VectorXd & getEeVelocityJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaV;
  }

  inline const Eigen::VectorXd & getEeVelocityJump()
  {
    //  std::cout<<"calling eevel jump: "<<getImpactBody()<<std::endl;
    const auto ee = cache_.grfContainer.find(getImpactBody());

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
  inline const Eigen::VectorXd & getBranchJointVelJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaQDot;
  }
  inline const Eigen::VectorXd & getBranchTauJump(const std::string & eeName) const
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.deltaTau;
  }

  inline void resetDataStructure()
  {
    // getOsd_()->resetDataStructure();
    cache_.reset();
  }
  inline void initializeDataStructure(int numEE)
  {

    cache_.ini(getOsd_()->getJacobianDim(), getOsd_()->getDof());
    // getOsd_()->initializeDataStructure(numEE);
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
  inline const Eigen::VectorXd & getImpulsiveForce()
  {
    const auto & ee = cache_.grfContainer.find(getImpactBody());
    return ee->second.impulseForce;
  }
  inline const sva::ForceVecd & getImpulsiveForceCOM()
  {
    /*
    sva::PTransformd X_ee_CoM = sva::PTransformd(getRobot().com())*getRobot().bodyPosW(getImpactBody()).inv();

    return X_ee_CoM.dualMul(sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce()));
    */
    const auto & ee = cache_.grfContainer.find(getImpactBody());
    return ee->second.impulseForceCOM;
  }
  const sva::ForceVecd & getImpulsiveForceCOP(const std::string & eeName);
  const sva::ForceVecd & getImpulsiveForceCOM(const std::string & eeName);
  const Eigen::VectorXd & getImpulsiveForce(const std::string & eeName);
  inline mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }
  inline void setImpactBody(std::string & impactBodyName)
  {
    impactBodyName_ = impactBodyName;
  }

  void setContact(std::string contactBodyName);
  // We need to specify the endeffectors that are in contact, e.g. two feet
  // std::vector<dart::dynamics::BodyNode *> contactEndEffectors;

  //----------------Public API for control ----------------------------/
  inline const Eigen::MatrixXd & getJacobianDeltaAlpha()
  {
    return cache_.jacobianDeltaAlpha;
  }
  inline const Eigen::MatrixXd & getJacobianDeltaTau()
  {
    return cache_.jacobianDeltaTau;
  }
  inline const Eigen::MatrixXd & getJacobianDeltaF(const std::string & eeName)
  {
    const auto & ee = cache_.grfContainer.find(eeName);
    return ee->second.jacobianDeltaF;
  }
  inline const double & getImpactDuration_() const
  {
    return impactDuration_;
  }

  inline const std::string & getImpactBody()
  {
    return impactBodyName_;
  }

protected:
  mc_rbdyn::Robot & robot_;
  const std::shared_ptr<mi_osd> & osdPtr_;

  std::string impactBodyName_;

  bool linearJacobian_;
  double impactDuration_;

  double coeFrictionDeduction_;
  double coeRes_;

  impactDataCache cache_;

  inline bool useLinearJacobian_() const
  {
    return linearJacobian_;
  }

  // impact end-effector
  // dart::dynamics::BodyNodePtr impactBodyPtr_;

  // The returned pointer is not supposed to be changed.
  // The returned pointer is not supposed to be changed.
  inline const std::shared_ptr<mi_osd> & getOsd_() const
  {
    return osdPtr_;
  }

  // mc_rbdyn::Robot & loadRobot_(mc_rbdyn::RobotModulePtr rm, const std::string & name);
  inline const double & getCoeRes_() const
  {
    return coeRes_;
  }
  inline const double & getCoeFricDe_() const
  {
    return coeFrictionDeduction_;
  }
  //  void tempTestEe_();
};
