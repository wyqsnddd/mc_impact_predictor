#pragma once

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotModule.h>

#include <RBDyn/FA.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/Jacobian.h>

#include <assert.h>
#include <map>

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
  mi_osd(mc_rbdyn::Robot & robot, bool linearJacobian);

  ~mi_osd() {}

  bool addEndeffector(std::string eeName);
  /*!
    \param EeNum the number of end-effectors to be used in the OSD.
    */
  void initializeDataStructure(int EeNum);

  void resetDataStructure();

  /*!
      \return the Operational space inertia matrix: \f$ \Lambda \f$
      */
  inline const Eigen::MatrixXd & getLambdaMatrix() const
  {
    return cache_.lambdaMatrix;
  }

  /*!
      \return the inverse of the Operational space inertia matrix: \f$ \Lambda^{-1} = JM^{-1}J^\top \f$
      */
  inline const Eigen::MatrixXd & getLambdaMatrixInv() const
  {
    return cache_.lambdaMatrixInv;
  }
  /*!
      \return the Operational space inertia matrix: \f$ \Lambda(i , j) \f$
      */
  const Eigen::MatrixXd getLambdaMatrix(int rowInt, int columnInt) const
  {
    return cache_.lambdaMatrix.block(rowInt * getEeNum(), columnInt * getEeNum(), 6, 6);
  }
  const int nameToIndex_(const std::string & eeName) const;

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
      \return
      \f$  \tilde{\Lambda}_{m}= ( \sum^{m}_{i=1}\Lambda_{mi}J_i ) \bar{J}_m \f$
  };
      */
  inline const Eigen::MatrixXd & getEffectiveLambdaMatrix(const std::string & eeName) const
  {
    return cache_.effectiveLambdaMatrices[nameToIndex_(eeName)];
  }
  /*
  Eigen::VectorXd getOsdForce(const std::string eeName) const
  {
    return cache_.osdF.segment(nameToIndex_(eeName) * getJacobianDim(), getJacobianDim());
  }
*/

  /*!
      \return  the dynamically consistent Jacobian inverse of the endeffector m : \f$ \bar{J}_m = ( (
     \sum^{m}_{i=1}\Lambda_{mi}J_i )M^{-1})^\top.\f$
      */
  inline const Eigen::MatrixXd & getDcJacobianInv(const std::string eeName) const
  {
    return cache_.dcJacobianInvs[nameToIndex_(eeName)];
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

  inline int getDof() const
  {
    return robotDof_;
  }
  inline int getEeNum() const
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
  inline const int getJacobianDim() const
  {
    return jacobianDim_;
  }

private:
  int robotDof_;
  bool linearJacobian_;
  const bool useLinearJacobian_()
  {
    return linearJacobian_;
  }
  int eeNum_;
  osdDataCache cache_;
  int jacobianDim_;
  mc_rbdyn::Robot & robot_;
  std::shared_ptr<rbd::ForwardDynamics> FDPtr_;

  void updateCache_();
};
