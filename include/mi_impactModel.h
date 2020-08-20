/* Copyright 2019 CNRS-UM LIRMM
 *
 * \author Yuquan Wang, Arnaud Tanguy 
 *
 * 
 *
 * mc_impact_predictor is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * mc_impact_predictor is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mc_impact_predictor. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <mc_rbdyn/Robot.h>

#include <RBDyn/Jacobian.h>

#include <Eigen/StdVector>
#include <iostream>

#include "mi_utils.h"

namespace mc_impact
{

class mi_impactModel
/*!
 * \brief Interface to the impact dynamics models.
 */
{
public:
  mi_impactModel(const mc_rbdyn::Robot & simRobot,
		 const ImpactModelParams & params)
  : simRobot_(simRobot), params_(params)
  {
    int dof = simRobot_.mb().nrDof();
    jacPtr_ = std::make_shared<rbd::Jacobian>(simRobot_.mb(), getParams().iBodyName );
    jacobian_.resize(getParams().dim, dof);
    jacobianDot_.resize(getParams().dim, dof);
    params_.inertial_surfaceNormal.normalize();
    robotJointVel_.resize(dof);
  }

  ~mi_impactModel() {}

  const ImpactModelParams & getParams() const
  {
    return params_; 
  }

  /*
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
  
  inline bool useBodyJacobian() const
  {
    return useBodyJacobian_; 
  }
  


*/
  inline const Eigen::VectorXd & getEeVelocity() const
  {
    return eeV_;
  }
  inline const Eigen::VectorXd & getEeVelocityJump() const
  {
    return deltaV_;
  }
  /*! 
   * @return P*(J + J_dot*dt)
   */
  inline const Eigen::MatrixXd & getProjectorTwo() const
  {
    return reductionProjectorTwo_;
  }

  /*! 
   * @return P*J
   */
  inline const Eigen::MatrixXd & getProjector() const
  {
    return reductionProjector_;
  }
  inline const Eigen::VectorXd & getJointVel() const
  {
    return robotJointVel_;
  }
  inline const Eigen::Vector3d & getSurfaceNormal()
  {
    return local_surfaceNormal_;
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

    inline const Eigen::MatrixXd & getJacobianDot()
  {
    return jacobianDot_;
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
  ImpactModelParams params_;
  // const std::shared_ptr<mi_osd> & osdPtr_;

  std::shared_ptr<rbd::Jacobian> jacPtr_;
  void updateJacobian_();
  Eigen::MatrixXd jacobian_;
  Eigen::MatrixXd jacobianDot_;


  ///< This is the impact normal direction in the inertial frame

  /*
  Eigen::Vector3d inertial_surfaceNormal_;
  std::string impactBodyName_;
  bool useBodyJacobian_ = true;
  double impactDuration_;
  double timeStep_;
  double coeFrictionDeduction_;
  double coeRes_;
  int dim_;
  */

  void update_();

  Eigen::VectorXd deltaV_ = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd eeV_ = Eigen::VectorXd::Zero(3);
  Eigen::MatrixXd reductionProjector_ = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd reductionProjectorTwo_ = Eigen::MatrixXd::Zero(3, 3);
  Eigen::VectorXd robotJointVel_;
  Eigen::Vector3d local_surfaceNormal_;
  Eigen::Vector3d contactVel_ = Eigen::Vector3d::Zero(3);
};
} // namespace mc_impact
