#pragma once

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotModule.h>

#include <RBDyn/FA.h>
#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/Jacobian.h>

//#include <dart/dynamics/dynamics.hpp>

#include <assert.h>
#include <map>


struct objJac{
	int containerIndex;  
	std::shared_ptr<rbd::Jacobian>	jacPtr;
};

struct osdDataCache
{
  // osdJacobian, osdJ_dot, dcJacobian, effectiveLambda, x_dot, tau have the same order.
  Eigen::MatrixXd osdJacobian;
  Eigen::MatrixXd osdJacobianDot;

  Eigen::MatrixXd invMassMatrix;
  Eigen::MatrixXd lambdaMatrix;
  // Eigen::MatrixXd crossLambdaMatrix;

  // Eigen::VectorXd osdAcc;
  //  Eigen::VectorXd osdVel;

  Eigen::MatrixXd lambdaMatrixInv;

  // Hash table: bodyNode <-> index in the local container
  //std::map<std::string, std::pair<std::shared_ptr<rbd::Jacobian>, int>> jacobians;
  std::map<std::string, objJac> jacobians;

  // This is a vector of the dynamically consistent Jacobian pseudo inverse of all the end-effectors.
  std::vector<Eigen::MatrixXd> dcJacobianInvs;
  // Eigen::MatrixXd dcJacobianInv;

  // This is a vector of the "effective operational mass matrix for each endeffector"
  std::vector<Eigen::MatrixXd> effectiveLambdaMatrices;
};

class mi_osd
{
public:
  mi_osd(mc_rbdyn::Robot & robot,
         // std::shared_ptr<rbd::ForwardDynamics> & fdPtr,
         bool linearJacobian);

  ~mi_osd() {}
  
  bool addEndeffector(std::string eeName);
  void initializeDataStructure(int EeNum);
  void resetDataStructure();

  Eigen::MatrixXd getLambdaMatrix() const
  {
    return cache_.lambdaMatrix;
  }
  Eigen::MatrixXd getLambdaMatrixInv() const
  {
    return cache_.lambdaMatrixInv;
  }
  const Eigen::MatrixXd getLambdaMatrix(int rowInt, int columnInt) const
  {
    return cache_.lambdaMatrix.block(rowInt * getEeNum(), columnInt * getEeNum(), 6, 6);
  }
  const int nameToIndex_(const std::string & eeName) const
  {
    auto tempEe = cache_.jacobians.find(eeName);
    if(tempEe != cache_.jacobians.end())
      //return tempEe->second.second;
      return tempEe->second.containerIndex;
    else
    {
      // std::cout << "Link " << eeName << " is missing." << std::endl;
      std::string error_msg = std::string("OSD::nameToIndex_: link-") + eeName + std::string(": does not exist.");
      throw std::runtime_error(error_msg);
    }
  }
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

  const Eigen::MatrixXd getJacobian(const std::string & eeName) const
  {
    return cache_.osdJacobian.block(nameToIndex_(eeName) * getJacobianDim(), 0, getJacobianDim(), getDof());
  }
  const Eigen::MatrixXd getJacobianDot(const std::string & eeName) const
  {
    return cache_.osdJacobianDot.block(nameToIndex_(eeName) * getJacobianDim(), 0, getJacobianDim(), getDof());
  }

  const Eigen::MatrixXd getEffectiveLambdaMatrix(const std::string & eeName) const
  {
    return cache_.effectiveLambdaMatrices[nameToIndex_(eeName)];
  }

  const Eigen::MatrixXd getDcJacobianInv(const std::string eeName) const
  {
    return cache_.dcJacobianInvs[nameToIndex_(eeName)];
  }
  const Eigen::MatrixXd getInvMassMatrix() const
  {
    return cache_.invMassMatrix;
  }
  // This needs to be called in every iteration only once
  void update()
  {
    if(getEeNum() < 3)
    {
      throw std::runtime_error("OSD: Too less end-effectors are defined for the OSD. ");
    }
    // mc_rtc components
    // std::cout << "Updating OSD FD..." << std::endl;

    // rbd::MultiBodyConfig & tempMbc = getRobot().mbc();

    rbd::forwardKinematics(getRobot().mb(), getRobot().mbc());
    rbd::forwardVelocity(getRobot().mb(), getRobot().mbc());
    rbd::forwardAcceleration(getRobot().mb(), getRobot().mbc());
    FDPtr_->forwardDynamics(getRobot().mb(), getRobot().mbc());
    // FDPtr_->computeH(getRobot().mb(), getRobot().mbc());
    // std::cout << "FD computed M ..." << std::endl;
    // std::cout << "Updating componentUpdateOsdDataCache_ ..." << std::endl;
    updateCache_();
  }
  int getDof() const
  {
    return robotDof_;
  }
  int getEeNum() const
  {
    return eeNum_;
  }
  mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }

  const std::shared_ptr<rbd::ForwardDynamics> getFD() const
  {
    return FDPtr_;
  }
  const int getJacobianDim() const
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
  // dart::dynamics::SkeletonPtr robotPtr_;
  mc_rbdyn::Robot & robot_;
  std::shared_ptr<rbd::ForwardDynamics> FDPtr_;

  /// Based on the symmetry, we calculate the inverse component wise.
  void updateCache_();
};
