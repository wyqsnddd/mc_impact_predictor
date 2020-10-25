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
#include <mc_rbdyn/RobotModule.h>

#include <RBDyn/FA.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/Jacobian.h>
#include <RBDyn/Momentum.h>

#include <assert.h>
#include <chrono>
#include <map>

namespace mc_impact
{

struct objJac
{
  int containerIndex;
  std::shared_ptr<rbd::Jacobian> jacPtr;
};

struct osdDataCache
/** \brief Data for computing the Operational Space equations of motion
 */
{
  // osdJacobian, osdJ_dot, dcJacobian, effectiveLambda, x_dot, tau have the same order.
  Eigen::MatrixXd osdJacobian; ///< Stack of Jacobians \f$ J = [J_1^T, J_2^T ....]^T \f$.
  Eigen::MatrixXd osdJacobianDot; ///< Derivative of osdJacobian

  Eigen::MatrixXd invMassMatrix; ///< Inverse of the joint space mass matrix
  Eigen::MatrixXd lambdaMatrix; ///< Operational space mass matrix

  // Eigen::VectorXd osdF;
  // Eigen::VectorXd osdAcc;
  //  Eigen::VectorXd osdVel;

  Eigen::MatrixXd lambdaMatrixInv; ///< Inverse of the OSD mass matrix

  // Hash table: bodyNode <-> index in the local container
  // std::map<std::string, std::pair<std::shared_ptr<rbd::Jacobian>, int>> jacobians;
  std::map<std::string, objJac> jacobians;

  // This is a vector of the dynamically consistent Jacobian pseudo inverse of all the end-effectors.
  std::vector<Eigen::MatrixXd>
      dcJacobianInvs; ///< the dynamically consistent Jacobian inverse: \f$ J_i = M^-1 J^T \Lambda \f$
  // Eigen::MatrixXd dcJacobianInv;

  // This is a vector of the "effective operational mass matrix for each endeffector"
  std::vector<Eigen::MatrixXd>
      effectiveLambdaMatrices; ///< \f$  \tilde{\Lambda}_{m}= ( \sum^{m}_{i=1}\Lambda_{mi}J_i ) \bar{J}_m \f$
};

class mi_osd
/*! \brief Operational Space Dynamics.
 */
{
public:
  /*!
      \param robot reference to the robot model used by the QP
      \param linearJacobian whether use the linear part of the Jacobian?
      \param bodyJacobian whether use the linear part of the Jacobian?
      */
  mi_osd(mc_rbdyn::Robot & robot, bool linearJacobian = true, bool bodyJacobian = true);

  ~mi_osd() {}

  void setContact(std::vector<std::string> & ees);

  bool hasEndeffector(const std::string & eeName);

  bool addEndeffector(std::string eeName);
  /*!
   * \param EeNum the number of end-effectors to be used in the OSD.
   */
  void initializeDataStructure(int EeNum);

  void resetDataStructure();

  /*!
   *  \return the operational-space inertia matrix: 
     \f$ \Lambda  = (J M^{-1} J^\top)^{-1}\f$
   */
  inline const Eigen::MatrixXd & getLambdaMatrix() const
  {
    return cache_.lambdaMatrix;
  }

  inline const std::vector<std::string> & getEes()
  {
    return endEffectors_;
  }

  inline const std::vector<std::string> & getContactEes()
  {
    return contactEndeffectors_;
  }
  inline std::size_t getContactNum()
  {
    return contactEndeffectors_.size();
  }
  /*!
   *  \return the inverse of the operational-space inertia matrix: 
   *  \f$ \Lambda^{-1} = JM^{-1}J^\top \f$
   */
  inline const Eigen::MatrixXd & getLambdaMatrixInv() const
  {
    return cache_.lambdaMatrixInv;
  }

  /*!
   * return the equivalent mass
   * \brief Robot arm (single kinematic chain): \f$ (J M^{-1} J^\top)^{-1} \f$
   * Humanoid: (kinematic tree): the corresponding block from the operational space inertia matrix: lambda.
   */
  inline const Eigen::MatrixXd getEquivalentMass(const std::string & eeName) const
  {
    return cache_.lambdaMatrix.block(nameToIndex_(eeName) * getJacobianDim(),
                                        nameToIndex_(eeName) * getJacobianDim(), getJacobianDim(), getJacobianDim());
  }

  /*!
    \return the Operational space inertia matrix: \f$ \Lambda(i , j) \f$
  */
  const Eigen::MatrixXd getLambdaMatrix(int rowInt, int columnInt) const
  {
    return cache_.lambdaMatrix.block(rowInt * static_cast<int>(getEeNum()), columnInt * static_cast<int>(getEeNum()), 6, 6);
  }
  const int & nameToIndex_(const std::string & eeName) const;

  const Eigen::MatrixXd getLambdaMatrix(const std::string & eeOne, const std::string & eeTwo) const
  {
    return cache_.lambdaMatrix.block(nameToIndex_(eeOne) * getJacobianDim(), nameToIndex_(eeTwo) * getJacobianDim(),
                                     getJacobianDim(), getJacobianDim());
  }
  const Eigen::MatrixXd getLambdaMatrixInv(const std::string & eeOne, const std::string & eeTwo) const
  {
    return cache_.lambdaMatrixInv.block(nameToIndex_(eeOne) * getJacobianDim(), nameToIndex_(eeTwo) * getJacobianDim(),
                                        getJacobianDim(), getJacobianDim());
  }

  const Eigen::MatrixXd getJacobian(const std::string & eeName)
  {
    return cache_.osdJacobian.block(nameToIndex_(eeName) * getJacobianDim(), 0, getJacobianDim(), getDof());
  }
  const Eigen::MatrixXd getJacobianDot(const std::string & eeName) const
  {
    return cache_.osdJacobianDot.block(nameToIndex_(eeName) * getJacobianDim(), 0, getJacobianDim(), getDof());
  }

  /*!
   * \return
   * \f$  \tilde{\Lambda}_{m}= ( \sum^{m}_{i=1}\Lambda_{mi}J_i ) \bar{J}_m \f$
   */
  inline const Eigen::MatrixXd & getEffectiveLambdaMatrix(const std::string & eeName) const
  {
    return cache_.effectiveLambdaMatrices[static_cast<unsigned long>(nameToIndex_(eeName))];
  }
  /*
  Eigen::VectorXd getOsdForce(const std::string eeName) const
  {
    return cache_.osdF.segment(nameToIndex_(eeName) * getJacobianDim(), getJacobianDim());
  }
*/

  /*!
      \return  the dynamically consistent Jacobian inverse of the endeffector m: 
      \f$ \bar{J}_m = ( (\sum^{m}_{i=1}\Lambda_{mi}J_i )M^{-1})^\top.\f$
      */
  inline const Eigen::MatrixXd & getDcJacobianInv(const std::string eeName) const
  {
    return cache_.dcJacobianInvs[static_cast<unsigned long>(nameToIndex_(eeName))];
  }

  inline const Eigen::MatrixXd & getMassMatrix() const
  {
    return getFD()->H();
  }

  /*!
      \return \f$ M^{-1}\f$
  */
  inline const Eigen::MatrixXd & getInvMassMatrix() const
  {
    return cache_.invMassMatrix;
  }

  /*!
   * This needs to be called in every iteration only once
   */
  void update();

  /*! \brief Time to update the internal dynamics models.
   * \return time in microseconds.
   */
  inline double modelUpdateTime() const
  {
    return modelUpdateTime_;
  }
  /*! \brief Time to solve the optimization problem.
   * \return time in microseconds.
   */
  inline double computationTime() const
  {
    return computationTime_;
  }
  inline int getDof() const
  {
    return robotDof_;
  }
  inline size_t getEeNum() const
  {
    return eeNum_;
  }
  inline mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }

  inline const std::shared_ptr<rbd::ForwardDynamics> getFD() const
  {
    return FDPtr_;
  }
  /*!
   * \return number of rows of the Jacobian
   */
  inline int getJacobianDim() const
  {
    return jacobianDim_;
  }
  void printInfo();
  bool useBodyJacobian() const
  {
    return bodyJacobian_;
  }
  bool useLinearJacobian() const
  {
    return linearJacobian_;
  }

  Eigen::Matrix3d crossMatrix(const Eigen::Vector3d & input);

  /*! \brief returns the 6 by 3 "force" (or linear) part of the Grasp matrix
   *  \param eeName of the end-effector where the force is applied
   *  \param reference point, e.g. the COM. By default, we use the origin of the inertial frame.
   *
   *  The reference frame is the inertial frame
   */
  Eigen::MatrixXd forceGraspMatrix(const std::string eeName,
                                   const Eigen::Vector3d & reference = Eigen::Vector3d::Zero());

  const sva::ForceVecd & getSimulatedCentroidalMomentum()
  {
    return centroidalMomentum_;
  }

  const sva::ForceVecd & getSimulatedCentroidalMomentumD()
  {
    return centroidalMomentumD_;
  }

private:
  int robotDof_;
  size_t eeNum_;
  osdDataCache cache_;
  int jacobianDim_;
  mc_rbdyn::Robot & robot_;

  bool linearJacobian_;
  bool bodyJacobian_;

  std::shared_ptr<rbd::ForwardDynamics> FDPtr_;

  std::vector<std::string> contactEndeffectors_; ///< end-effectors with established contact.
  std::vector<std::string> endEffectors_; ///< end-effectors

  bool addEndeffector_(std::string eeName);
  void updateCache_();

  double computationTime_;
  double modelUpdateTime_;

  sva::ForceVecd centroidalMomentumD_ = sva::ForceVecd::Zero();
  sva::ForceVecd centroidalMomentum_ = sva::ForceVecd::Zero();
};
} // namespace mc_impact
