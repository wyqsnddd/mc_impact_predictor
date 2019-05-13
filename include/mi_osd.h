#pragma once

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotModule.h>

#include <RBDyn/FD.h>
#include <RBDyn/FA.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/Jacobian.h>

//#include <dart/dynamics/dynamics.hpp>

#include <assert.h>
#include <map>

struct osdDataCache
{
  // osdJacobian, osdJ_dot, dcJacobian, effectiveLambda, x_dot, tau have the same order.
  Eigen::MatrixXd osdJacobian;
  Eigen::MatrixXd osdJacobianDot;

  Eigen::MatrixXd invMassMatrix;
  Eigen::MatrixXd lambdaMatrix;

  Eigen::VectorXd rhoOne;
  Eigen::VectorXd rhoTwo;

  Eigen::VectorXd osdAcc;
  Eigen::VectorXd osdVel;
  Eigen::VectorXd osdTau;

  Eigen::MatrixXd lambdaMatrixInv;


  // Hash table: bodyNode <-> index in the local container
  std::map<std::string, std::pair<std::shared_ptr<rbd::Jacobian>, int>> jacobians;

  // This is a vector of the dynamically consistent Jacobian pseudo inverse of all the end-effectors.
  std::vector<Eigen::MatrixXd> dcJacobianInvs;

  // This is a vector of the "effective operational mass matrix for each endeffector"
  std::vector<Eigen::MatrixXd> effectiveLambdaMatrices;
};

class mi_osd
{
public:
  mi_osd(// const dart::dynamics::SkeletonPtr & robotPtr,
		  mc_rbdyn::Robot & robot, 
		  bool linearJacobian);

  ~mi_osd() {}

  Eigen::MatrixXd getLambdaMatrix() const
  {
    return cache_.lambdaMatrix;
  }
  Eigen::MatrixXd getLambdaMatrixInv() const
  {
    return cache_.lambdaMatrixInv;
  }
  Eigen::MatrixXd getLambdaMatrix(int rowInt, int columnInt) const
  {
    return cache_.lambdaMatrix.block(rowInt * getDof(), columnInt * getDof(), 6, 6);
  }
  const Eigen::MatrixXd getJacobian(const std::string & eeName) const
  {
    std::cout << "Osd looks for Jacobian of " << eeName << std::endl;
    /*
    for(auto it = cache_.jacobians.begin(); it != cache_.jacobians.end(); ++it)
    {

      std::cout << it->first << " has a local index: " << it->second.second << std::endl<<
    cache_.osdJacobian.block(it->second.second* jacobianDim_, 0, jacobianDim_, getDof())<<std::endl;
      //int ii = it->second.second;
    }
    */
    /*
    if(cache_.jacobians.find(eeName) != cache_.jacobians.end())
    {
      std::cout << "Key found";
    }
    else
    {
      std::cout << "Key not found";
    }
    */
    auto ee = cache_.jacobians.find(eeName);

    std::cout << "Osd found ee" << ee->first << std::endl;

    int index = ee->second.second;
    std::cout << "Endeffector: " << ee->first << " has index: " << index << std::endl;
    return cache_.osdJacobian.block(index * jacobianDim_, 0, jacobianDim_, getDof());
  }
  const Eigen::MatrixXd getEffectiveLambdaMatrix(const std::string & eeName) const
  {
    auto ee = cache_.jacobians.find(eeName);
    return cache_.effectiveLambdaMatrices[ee->second.second];
  }

  const Eigen::MatrixXd getDcJacobianInv(const std::string eeName) const
  {
    auto ee = cache_.jacobians.find(eeName);
    return cache_.dcJacobianInvs[ee->second.second];
  }
  const Eigen::MatrixXd getInvMassMatrix() const
  {
    return cache_.invMassMatrix;
  }
  // This needs to be called in every iteration only once
  void update()
  {
    // mc_rtc components
    std::cout << "Updating OSD FD..." << std::endl;
    //FDPtr_->forwardDynamics(getRobot().mb(), const_cast<rbd::MultiBodyConfig & >(getRobot().mbc()));
    
    //rbd::MultiBodyConfig & tempMbc = getRobot().mbc();

    rbd::forwardKinematics(getRobot().mb(), getRobot().mbc() );
    rbd::forwardVelocity(getRobot().mb(), getRobot().mbc() );
    rbd::forwardAcceleration(getRobot().mb(), getRobot().mbc());
    FDPtr_->forwardDynamics(getRobot().mb(), getRobot().mbc() );
    //FDPtr_->computeH(getRobot().mb(), getRobot().mbc());
    std::cout << "FD computed M ..." << std::endl;
    std::cout << "Updating componentUpdateOsdDataCache_ ..." << std::endl;
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
  /*
  const mc_rbdyn::Robot & getRobot() 
  {
    return robot_;
  }
  */
  mc_rbdyn::Robot & getRobot() 
  {
    return robot_;
  }
  const std::shared_ptr<rbd::ForwardDynamics> getFD() const
  {
    return FDPtr_;
  }
  bool nonSingular()
  {
    return nonSingular_;
  }

private:
  int robotDof_;
  bool linearJacobian_;
  bool useLinearJacobian_()
  {
    return linearJacobian_;
  }
  int eeNum_;
  osdDataCache cache_;
  int jacobianDim_;
  bool nonSingular_;
/*
  const dart::dynamics::SkeletonPtr getDartRobot() const{
    return robotPtr_; 
  }
  std::size_t getDartBodyIndex_(const std::string input) const{
    return getDartRobot()->getIndexOf(getDartRobot()->getBodyNode(input));
  }
*/
  //dart::dynamics::SkeletonPtr robotPtr_;
  mc_rbdyn::Robot & robot_;
  std::shared_ptr<rbd::ForwardDynamics> FDPtr_;
  /// Direct inverse
  //void intuitiveUpdateOsdDataCache_();

  /// Based on the symmetry, we calculate the inverse component wise.
  void updateCache_();
};
