#pragma once

#include <mc_rbdyn/Robot.h>

#include <RBDyn/Jacobian.h>

#include <Eigen/StdVector>
#include <iostream>

namespace mc_impact
{

class mi_impactModel
{
public:
  mi_impactModel(const mc_rbdyn::Robot & simRobot,
                 // const std::shared_ptr<mi_osd> & osdPtr,
                 const std::string & iBodyName,
                 const Eigen::Vector3d & inertial_surfaceNormal,
                 double iDuration = 0.005,
                 double timeStep = 0.005,
                 double coeF = 0.2,
                 double coeR = 0.8,
                 int dim = 3)
  : simRobot_(simRobot), impactBodyName_(iBodyName), inertial_surfaceNormal_(inertial_surfaceNormal),
    impactDuration_(iDuration), timeStep_(timeStep), coeFrictionDeduction_(coeF), coeRes_(coeR), dim_(dim)
  {
    jacPtr_ = std::make_shared<rbd::Jacobian>(simRobot_.mb(), getImpactBody());
    jacobian_.resize(getDim(), simRobot_.mb().nrDof());
    inertial_surfaceNormal_.normalize();
  }

  ~mi_impactModel() {}

  inline int getDim() const
  {
    return dim_;
  }
  inline const std::string & getImpactBody() const
  {
    return impactBodyName_;
  }
  inline double getTimeStep() const
  {
    return timeStep_;
  }
  inline double getImpactDuration() const
  {
    return impactDuration_;
  }
  inline double getCoeRes() const
  {
    return coeRes_;
  }
  inline double getCoeFricDe() const
  {
    return coeFrictionDeduction_;
  }
  inline const Eigen::VectorXd & getEeVelocity() const
  {
    return eeV_;
  }
  inline const Eigen::VectorXd & getEeVelocityJump() const
  {
    return deltaV_;
  }
  inline const Eigen::MatrixXd & getProjector() const
  {
    return reductionProjector_;
  }
  inline const Eigen::VectorXd & getJointVel() const
  {
    return temp_q_vel_;
  }

  /** Predict the impact-induced state jumps based on the internal update.
   *
   */
  void update();

  /** Predict the impact-induced state jumps based on given impact normal direction
   *
   * @param surfaceNormal the given impact normal direction
   */
  void update(const Eigen::Vector3d & surfaceNormal);
  inline const Eigen::Vector3d & getSurfaceNormal()
  {
    return local_surfaceNormal_;
  }

  inline const Eigen::MatrixXd & getJacobian()
  {
    return jacobian_;
  }
  inline const Eigen::Vector3d & getContactVel()
  {
    return contactVel_; 
  }
private:
  const mc_rbdyn::Robot & simRobot_;
  // const std::shared_ptr<mi_osd> & osdPtr_;

  std::shared_ptr<rbd::Jacobian> jacPtr_;
  void updateJacobian_();
  Eigen::MatrixXd jacobian_;

  std::string impactBodyName_;

  ///< This is the impact normal direction in the inertial frame
  Eigen::Vector3d inertial_surfaceNormal_;

  double impactDuration_;
  double timeStep_;
  double coeFrictionDeduction_;
  double coeRes_;
  int dim_;

  void update_();

  Eigen::VectorXd deltaV_ = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd eeV_ = Eigen::VectorXd::Zero(3);
  Eigen::MatrixXd reductionProjector_ = Eigen::MatrixXd::Zero(3, 3);
  Eigen::VectorXd temp_q_vel_;
  Eigen::Vector3d local_surfaceNormal_;
  Eigen::Vector3d contactVel_ = Eigen::Vector3d::Zero(3);
};
} // namespace mc_impact
