#include "mi_osd.h"

mi_osd::mi_osd(mc_rbdyn::Robot & robot,
               //    std::shared_ptr<rbd::ForwardDynamics> & fdPtr,
               bool linearJacobian)
: robot_(robot) //, FDPtr_(fdPtr) // robotPtr_(robotPtr),
{
  std::cout << "The osd dynamics constructor is called " << std::endl;
  linearJacobian_ = linearJacobian;
  // Initilize the forward dynamics:
  FDPtr_ = std::make_shared<rbd::ForwardDynamics>(getRobot().mb());

  std::cout << "The FD constructor is built." << std::endl;
  // Create a local copy to avoid touching the mc_rtc controller robot.
  // rbd::MultiBodyConfig & tempMbc = getRobot().mbc();
  rbd::forwardKinematics(getRobot().mb(), getRobot().mbc());
  rbd::forwardVelocity(getRobot().mb(), getRobot().mbc());
  rbd::forwardAcceleration(getRobot().mb(), getRobot().mbc());
  FDPtr_->forwardDynamics(getRobot().mb(), getRobot().mbc());
  // FDPtr_->computeH(getRobot().mb(), getRobot().mbc());
  // Initialize the Jacobians
  int mRows = static_cast<int>(getFD()->H().rows());
  int mCols = static_cast<int>(getFD()->H().cols());
  assert(mCols == mRows);
  robotDof_ = mRows;

  Eigen::MatrixXd tempMassMatrix = getFD()->H();
  // std::cout<<"The mass matrix is: " <<getFD()->H()<<std::endl;

  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_mM(tempMassMatrix);
  cache_.invMassMatrix.resize(mRows, mCols);
  cache_.invMassMatrix = lu_decomp_mM.inverse();

  // assert(mRows == mCols);

  if(useLinearJacobian_())
  {
    jacobianDim_ = 3;
  }
  else
  {
    jacobianDim_ = 6;
  }
  // getJacobianDim() = 6;

  // std::cout << "Updated OSD." << std::endl;
  cache_.jacobians = std::map<std::string, std::pair<std::shared_ptr<rbd::Jacobian>, int>>();
  eeNum_ = 0;
  std::cout << "OSD is created." << std::endl;
}
void mi_osd::resetDataStructure()
{
  cache_.osdJacobian.setZero();
  cache_.osdJacobianDot.setZero();
  // cache_.osdAcc.setZero();
  // cache_.osdVel.setZero();
  cache_.lambdaMatrix.setZero();
  cache_.lambdaMatrixInv.setZero();
  cache_.invMassMatrix.setZero();
  for(int ii = 0; ii < getEeNum(); ii++)
  {
    cache_.dcJacobianInvs[ii].resize(getDof(), getJacobianDim());
    cache_.effectiveLambdaMatrices[ii].resize(getJacobianDim(), getDof());
  }
}
void mi_osd::initializeDataStructure(const int EeNum)
{
  // eeNum_ = static_cast<int>(cache_.jacobians.size());

  cache_.osdJacobian.resize(EeNum * getJacobianDim(), getDof());
  cache_.osdJacobianDot.resize(EeNum * getJacobianDim(), getDof());
  // cache_.osdAcc.resize(EeNum * getJacobianDim());
  // cache_.osdVel.resize(EeNum * getJacobianDim());
  cache_.lambdaMatrix.resize(EeNum * getJacobianDim(), EeNum * getJacobianDim());
  cache_.lambdaMatrixInv.resize(EeNum * getJacobianDim(), EeNum * getJacobianDim());

  cache_.dcJacobianInvs.resize(EeNum);
  cache_.effectiveLambdaMatrices.resize(EeNum);
  /// Get the robot end-effectors
  for(int ii = 0; ii < getEeNum(); ii++)
  {
    cache_.dcJacobianInvs[ii].resize(getDof(), getJacobianDim());
    cache_.effectiveLambdaMatrices[ii].resize(getJacobianDim(), getDof());
  }
}
void mi_osd::updateCache_()
{
  // Read from the robot:
  //  std::cout << "Updating OSD cache..." << std::endl;
  // Eigen::MatrixXd tempMassMatrix = getFD()->H();
  // std::cout<<"The mass matrix is: " <<getFD()->H()<<std::endl;

  // std::cout<<"The C is: "<<getFD()->C()<<std::endl;
  // (0) Update the mass matrix inverse
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_M(getFD()->H());
  cache_.invMassMatrix = lu_decomp_M.inverse();

  // (1) Update the OSD Jacobian
  for(auto it = cache_.jacobians.begin(); it != cache_.jacobians.end(); ++it)
  {
    int ii = it->second.second;
    Eigen::MatrixXd tempJacobian = it->second.first->jacobian(getRobot().mb(), getRobot().mbc());
    Eigen::MatrixXd tempJacobianDot = it->second.first->jacobianDot(getRobot().mb(), getRobot().mbc());

    Eigen::MatrixXd tempFullJacobian, tempFullJacobianDot;
    tempFullJacobian.resize(getJacobianDim(), getDof());
    tempFullJacobianDot.resize(getJacobianDim(), getDof());

    if(useLinearJacobian_())
    {
      it->second.first->fullJacobian(getRobot().mb(), tempJacobian.block(3, 0, 3, tempJacobian.cols()),
                                     tempFullJacobian);
      it->second.first->fullJacobian(getRobot().mb(), tempJacobian.block(3, 0, 3, tempJacobianDot.cols()),
                                     tempFullJacobianDot);
    }
    else
    {
      it->second.first->fullJacobian(getRobot().mb(), tempJacobian, tempFullJacobian);
      it->second.first->fullJacobian(getRobot().mb(), tempJacobianDot, tempFullJacobianDot);
    }

    cache_.osdJacobian.block(ii * getJacobianDim(), 0, getJacobianDim(), getDof()) = tempFullJacobian;
    cache_.osdJacobianDot.block(ii * getJacobianDim(), 0, getJacobianDim(), getDof()) = tempFullJacobianDot;
  }

  // (2) Update the OSD Innertia matrix
  cache_.lambdaMatrixInv = cache_.osdJacobian * getInvMassMatrix() * (cache_.osdJacobian.transpose());

  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_lambda_inv(cache_.lambdaMatrixInv);
  cache_.lambdaMatrix = lu_decomp_lambda_inv.inverse();

  // auto tempJointTorque = rbd::dofToVector(getRobot().mb(), getRobot().mbc().jointTorque);

  // (3) Update the dynamically consistent Jacobian inverse:
  for(int ii = 0; ii < getEeNum(); ii++)
  {
    cache_.effectiveLambdaMatrices[ii] =
        cache_.lambdaMatrix.block(ii * getJacobianDim(), 0, getJacobianDim(), getEeNum() * getJacobianDim())
        * cache_.osdJacobian;

    cache_.dcJacobianInvs[ii] = (cache_.effectiveLambdaMatrices[ii] * getInvMassMatrix()).transpose();
  }
}
