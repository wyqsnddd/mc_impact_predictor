#include "mi_osd.h"

mi_osd::mi_osd(const dart::dynamics::SkeletonPtr & robotPtr,
               // const mc_rbdyn::Robot & robot,
               bool linearJacobian)
: // robot_(robot)
  robotPtr_(robotPtr)
{
  std::cout << "The osd dynamics constructor is called " << std::endl;
  linearJacobian_ = linearJacobian;
  // Initilize the forward dynamics:
  // FDPtr_ = std::make_shared<rbd::ForwardDynamics>(getRobot().mb());

  // std::cout << "The FD constructor is built." << std::endl;
  // Create a local copy to avoid touching the mc_rtc controller robot.
  // rbd::MultiBodyConfig tempMbc = getRobot().mbc();
  // FDPtr_->forwardDynamics(getRobot().mb(), tempMbc);
  // FDPtr_->computeH(getRobot().mb(), getRobot().mbc());
  // std::cout << "The masss matrix is built." << std::endl;
  // Initialize the endEffectors
  int mRows = static_cast<int>(getDartRobot()->getMassMatrix().rows());
  int mCols = static_cast<int>(getDartRobot()->getMassMatrix().cols());
  assert(mRows == mCols);

  /*
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_M(getFD()->H());
  cache_.invMassMatrix.resize(mRows, mCols);
  cache_.invMassMatrix = lu_decomp_M.inverse();
  */
  assert(mRows == mCols);
  robotDof_ = mRows;

  if(useLinearJacobian_())
  {
    jacobianDim_ = 3;
  }
  else
  {
    jacobianDim_ = 6;
  }
  // jacobianDim_ = 6;
  nonSingular_ = true;

  std::cout << "The stack of endEffectors are about to be built." << std::endl;
  // I can use a configuration file to specify the end-effectors later.
  /*
  cache_.endEffectors.insert(
      std::make_pair(
        strdup("l_ankle"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "l_ankle" ),
          0
          )
        )
      );

  cache_.endEffectors.insert(
      std::make_pair(
        strdup("r_ankle"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "r_ankle" ),
          1
          )
        )
      );

  cache_.endEffectors.insert(
      std::make_pair(
        strdup("l_wrist"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "l_wrist" ),
          2
          )
        )
      );
  cache_.endEffectors.insert(
      std::make_pair(
        strdup("r_wrist"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "r_wrist" ),
          3
          )
        )
      );
*/

  cache_.endEffectors["l_ankle"] = std::make_pair(getDartRobot()->getBodyNode("l_ankle"), 0);
  cache_.endEffectors["r_ankle"] = std::make_pair(getDartRobot()->getBodyNode("r_ankle"), 1);
  cache_.endEffectors["l_wrist"] = std::make_pair(getDartRobot()->getBodyNode("l_wrist"), 2);
  cache_.endEffectors["r_wrist"] = std::make_pair(getDartRobot()->getBodyNode("r_wrist"), 3);

  std::cout << "The stack of endEffectors are built." << std::endl;

  eeNum_ = static_cast<int>(cache_.endEffectors.size());

  cache_.osdJacobian.resize(getEeNum() * jacobianDim_, getDof());
  cache_.osdJacobianDot.resize(getEeNum() * jacobianDim_, getDof());
  cache_.osdAcc.resize(getEeNum() * jacobianDim_);
  cache_.osdVel.resize(getEeNum() * jacobianDim_);
  cache_.osdTau.resize(getEeNum() * jacobianDim_);
  /*
    contactEndEffectors.insert(std::make_pair(getRobot()->getBodyNode("l_ankle"), true));
    contactEndEffectors.insert(std::make_pair(getRobot()->getBodyNode("r_ankle"), true));
    contactEndEffectors.insert(std::make_pair(getRobot()->getBodyNode("l_wrist"), false));
    contactEndEffectors.insert(std::make_pair(getRobot()->getBodyNode("r_wrist"), false));
  */
  /// Initialize the cache:
  cache_.lambdaMatrix.resize(getEeNum() * jacobianDim_, getEeNum() * jacobianDim_);
  cache_.lambdaMatrixInv.resize(getEeNum() * jacobianDim_, getEeNum() * jacobianDim_);

  cache_.rhoOne.resize(getEeNum() * jacobianDim_);
  cache_.rhoTwo.resize(getEeNum() * jacobianDim_);

  cache_.dcJacobianInvs.resize(getEeNum());
  cache_.effectiveLambdaMatrices.resize(getEeNum());
  /// Get the robot end-effectors
  for(int ii = 0; ii < getEeNum(); ii++)
  {
    cache_.dcJacobianInvs[ii].resize(getDof(), jacobianDim_);
    cache_.effectiveLambdaMatrices[ii].resize(jacobianDim_, getDof());
  }

  std::cout << "Updating OSD..." << std::endl;
  update();

  std::cout << "Updated OSD." << std::endl;
}

void mi_osd::createFullJacobian(const Eigen::MatrixXd & inputJac,
                                Eigen::MatrixXd & fullJac,
                                const dart::dynamics::BodyNodePtr & bnPtr)
{

  fullJac.setZero();
  auto dependentDofs = bnPtr->getDependentGenCoordIndices();
  for(std::size_t ii = 0; ii < static_cast<size_t>(dependentDofs.size()); ii++)
  {

    std::size_t index = dependentDofs[ii];
    assert(bnPtr->dependsOn(index));
    fullJac.block(0, index, getJacobianDim_(), 1) = inputJac.block(0, ii, getJacobianDim_(), 1);
  }
}

void mi_osd::updateCache_()
{
  // Read from the robot:
  std::cout << "Updating OSD cache..." << std::endl;

  // FDPtr_->forwardDynamics(getRobot().mb(), getRobot().mbc());
  // FDPtr_->computeH(getRobot().mb(), getRobot().mbc());
  // Update the mass matrix inverse
  // Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_M(getFD()->H());
  // cache_.invMassMatrix = lu_decomp_M.inverse();
  std::cout << "Inverse mass matrix calculated, there are row: " << getDartRobot()->getInvMassMatrix().rows()
            << ", col: " << getDartRobot()->getInvMassMatrix().cols() << std::endl;

  for(auto it = cache_.endEffectors.begin(); it != cache_.endEffectors.end(); ++it)
  {
    int ii = it->second.second;
    std::cout << it->first << " has a local index: " << it->second.second << std::endl;
    std::cout << it->first << " has number of dependent joints: " << it->second.first->getNumDependentGenCoords()
              << std::endl;
    // cache_.osdJacobian.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) = it->first->getJacobian();
    //
    // auto tempJacobian = it->second.first->jacobian(getRobot().mb(), getRobot().mbc());
    // auto tempJacobianDot = it->second.first->jacobianDot(getRobot().mb(), getRobot().mbc());
    // std::cout<<"tempJacobian size is: "<<tempJacobian.rows()<<", "<<tempJacobian.cols()<<std::endl;
    /*
    std::cout << "Jacobian matrix calculated..." << std::endl;
    Eigen::MatrixXd tempFullJacobian, tempFullJacobianDot;
    tempFullJacobian.resize(jacobianDim_, getDof());
    tempFullJacobianDot.resize(jacobianDim_, getDof());
*/
    // std::cout<<"temp Full Jacobian size is: "<<tempFullJacobian.rows()<<", "<<tempFullJacobian.cols()<<std::endl;
    if(useLinearJacobian_())
    {
      /*
      it->second.first->fullJacobian(getRobot().mb(),
          tempJacobian.block(3, 0, 3, tempJacobian.cols()),
                                     tempFullJacobian);
      it->second.first->fullJacobian(getRobot().mb(),
          tempJacobian.block(3, 0, 3, tempJacobianDot.cols()),
                                     tempFullJacobianDot);
             */
      /*
      cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) =
        getDartRobot()->getBodyNode(it->first)->getCOMLinearAcceleration();
        */

      // std::cout<<"test 1"<<std::endl;
      /*
     auto jac_one = it->second.first->getWorldJacobian(Eigen::Vector3d::Zero()) ;
     auto jac_two = it->second.first->getJacobian(Eigen::Vector3d::Zero()) ;
     auto jac_three = it->second.first->getLinearJacobian(Eigen::Vector3d::Zero()) ;

      std::cout<<"The world jac has rows: "<<jac_one.rows()<<", and cols: "<<jac_one.cols()<<", it should have cols:
     "<<getDof()<<std::endl; std::cout<<"The jacobian jac has rows: "<<jac_two.rows()<<", and cols:
     "<<jac_two.cols()<<", it should have cols: "<<getDof()<<std::endl; std::cout<<"The linear jac has rows:
     "<<jac_three.rows()<<", and cols: "<<jac_three.cols()<<", it should have cols: "<<getDof()<<std::endl;

      //auto tempJac = it->second.first->(Eigen::Vector3d::Zero()).block(getJacobianDim_(), 0, getJacobianDim_(),
     getDof());
*/
      auto tempJac = it->second.first->getLinearJacobian(Eigen::Vector3d::Zero());

      Eigen::MatrixXd linearJacOne;
      linearJacOne.resize(getJacobianDim_(), getDof());
      createFullJacobian(tempJac, linearJacOne, it->second.first);
      // getLinearJacobian(dart::dynamics::Frame::World());
      //	    std::cout<<"The linear jac has rows: "<<linearJac.rows()<<", and cols: "<<linearJac.cols()<<", it should
      //have cols: "<<getDof()<<std::endl;
      cache_.osdJacobian.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) = linearJacOne;

      auto tempJacDot = it->second.first->getLinearJacobianDeriv(dart::dynamics::Frame::World());

      Eigen::MatrixXd linearJacDot;
      linearJacDot.resize(getJacobianDim_(), getDof());
      createFullJacobian(tempJacDot, linearJacDot, it->second.first);

      //	    std::cout<<"test 2"<<std::endl;
      cache_.osdJacobianDot.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) = linearJacDot;

      //	    std::cout<<"test 3"<<std::endl;
      cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) =
          getDartRobot()->getBodyNode(it->first)->getLinearAcceleration();

      //	    std::cout<<"test 4"<<std::endl;
      cache_.osdVel.segment(ii * jacobianDim_, jacobianDim_) =
          getDartRobot()->getBodyNode(it->first)->getLinearVelocity();
    }
    else
    {
      /*
      it->second.first->fullJacobian(getRobot().mb(), tempJacobian, tempFullJacobian);
      it->second.first->fullJacobian(getRobot().mb(), tempJacobianDot, tempFullJacobianDot);
      */
      /*
            cache_.osdJacobian.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) =
         it->second.first->getWorldJacobian(Eigen::Vector3d::Zero());

            cache_.osdJacobianDot.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) =
         it->second.first->getJacobianSpatialDeriv(dart::dynamics::Frame::World());

            cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) =
              getDartRobot()->getBodyNode(it->first)->getSpatialAcceleration();

            cache_.osdVel.segment(ii * jacobianDim_, jacobianDim_) =
         getDartRobot()->getBodyNode(it->first)->getSpatialVelocity();
            */
      /*
            cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) =
            getDartRobot()->getBodyNode(it->first)->getCOMSpatialAcceleration();

      */
      /*
      cache_.osdVel.segment(ii * jacobianDim_, jacobianDim_) =
            getRobot().mbc().bodyVelW
            [
            getRobot().mb().bodyIndexByName(it->first)
            ].vector();
*/
    }

    // std::cout << "Full Jacobian matrix calculated..." << std::endl;

    // std::cout << "" << std::endl;
    /*
 cache_.osdTau.block(ii * jacobianDim_, 0, jacobianDim_, 1) =
   getDartRobot()->getForce(
       getDartBodyIndex_(it->first)
       );
       */

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_J_temp(
        cache_.osdJacobian.block(ii * getJacobianDim_(), 0, getJacobianDim_(), getDof()));
    std::cout << it->first << " Jacobian has rank: " << lu_decomp_J_temp.rank() << std::endl;
  }

  std::cout << "OSD Jacobian calculated..." << std::endl;

  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_J(cache_.osdJacobian);
  std::cout << "The stacked OSD Jacobian has row: " << cache_.osdJacobian.rows()
            << ", col: " << cache_.osdJacobian.cols() << ", rank: " << lu_decomp_J.rank() << std::endl;

  cache_.lambdaMatrixInv = cache_.osdJacobian * getDartRobot()->getInvMassMatrix() * (cache_.osdJacobian.transpose());

  std::cout << "Inverse Lambda mass matrix calculated..." << std::endl;
  // Compute the Lambda matrix component-wise
  /*
  for(int ii = 0; ii < getEeNum(); ii++)
  {
    for(int jj = 0; jj < getEeNum(); jj++)
    {
      // std::cout<<"calculating row: "<<ii<<", col: "<<jj<<std::endl;
      cache_.lambdaMatrix.block(ii * jacobianDim_, jj * jacobianDim_, jacobianDim_, jacobianDim_) =
          (cache_.osdJacobian.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) * getInvMassMatrix()
           * cache_.osdJacobian.block(jj * jacobianDim_, 0, jacobianDim_, getDof()).transpose())
              .inverse();
    }
  }
  */

  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_lambda_inv(cache_.lambdaMatrixInv);
  cache_.lambdaMatrix = lu_decomp_lambda_inv.inverse();

  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_lambda(cache_.lambdaMatrix);
  std::cout << "The OSD mass matrix has row: " << cache_.lambdaMatrix.rows() << ", col: " << cache_.lambdaMatrix.cols()
            << ", rank: " << lu_decomp_lambda.rank() << std::endl;
  // std::cout<<"The OSD MassMatrix is: "<< std::endl<<cache_.lambdaMatrix<<std::endl;

  std::cout << "The joint torque is: " << std::endl << getDartRobot()->getForces() << std::endl;
  std::cout << "The number of ee is: " << getEeNum() << std::endl;

  // Update the dynamically consistent Jacobian inverse:
  for(int ii = 0; ii < getEeNum(); ii++)
  {
    // std::cout<<"ii: "<<ii<<std::endl;
    cache_.effectiveLambdaMatrices[ii] =
        cache_.lambdaMatrix.block(ii * getJacobianDim_(), 0, getJacobianDim_(), getEeNum() * getJacobianDim_())
        * cache_.osdJacobian;

    cache_.dcJacobianInvs[ii] = (cache_.effectiveLambdaMatrices[ii] * getInvMassMatrix()).transpose();

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_dcJ(cache_.dcJacobianInvs[ii]);

    // Update the OSD components:
    cache_.rhoTwo.segment(ii * getJacobianDim_(), getJacobianDim_()) =
        cache_.dcJacobianInvs[ii].transpose() * getDartRobot()->getCoriolisAndGravityForces();
    std::cout<<"The left ankle force is: "<< getDartRobot()->getBodyNode(getDartBodyIndex_("l_ankle"))->getBodyForce() <<std::endl;

    if(useLinearJacobian_())
    {
      cache_.osdTau.segment(ii * jacobianDim_, jacobianDim_) =
          cache_.dcJacobianInvs[ii].transpose()
          * (getDartRobot()->getForces()
             + getJacobian("l_ankle").transpose()
                   * getDartRobot()->getBodyNode(getDartBodyIndex_("l_ankle"))->getBodyForce().tail<3>()
             + getJacobian("r_ankle").transpose()
                   * getDartRobot()->getBodyNode(getDartBodyIndex_("r_ankle"))->getBodyForce().tail<3>());
    }
    else
    {
      cache_.osdTau.segment(ii * jacobianDim_, jacobianDim_) =
          cache_.dcJacobianInvs[ii].transpose()
          * (getDartRobot()->getForces()
             + getJacobian("l_ankle").transpose()
                   * getDartRobot()->getBodyNode(getDartBodyIndex_("l_ankle"))->getBodyForce()
             + getJacobian("r_ankle").transpose()
                   * getDartRobot()->getBodyNode(getDartBodyIndex_("r_ankle"))->getBodyForce());
    }
    std::cout << "The dynamically consistent Jacobian " << ii << " has row: " << cache_.dcJacobianInvs[ii].rows()
              << ", col: " << cache_.dcJacobianInvs[ii].cols() << ", rank: " << lu_decomp_dcJ.rank() << std::endl;
  }
  // Update the rest of the OSD components

  std::cout << "Rho One is to be calculated.";

  cache_.rhoOne = -cache_.lambdaMatrix * cache_.osdJacobianDot * getDartRobot()->getVelocities();

  std::cout << "Rho One is calculated.";

  /*
    cache_.rhoTwo = cache_.dcJacobianInvs.transpose()
        *getFD()->C();

    cache_.osdF = cache_.dcJacobianInvs.transpose()
        *cache_.osdTau;

  */

  // (0) This is a temporary test:

  std::cout << "Inverse of the OSD inertia matrix is: " << std::endl << cache_.lambdaMatrixInv << std::endl;
  std::cout << "The OSD inertia matrix is: " << std::endl << cache_.lambdaMatrix << std::endl;
  std::cout << "Multipulication of Lambda*Lambda_inv is: " << std::endl
            << cache_.lambdaMatrix * cache_.lambdaMatrixInv << std::endl;
  // (1) dc Jacobian inverse multiplies J should be close to identity?

  for(auto it = cache_.endEffectors.begin(); it != cache_.endEffectors.end(); ++it)
  {

    int index = it->second.second;
    // J*dc_J_inv:
    std::cout << "Test " << it->first << " Jacobian multiplies DC Jacobian inverse: " << std::endl
              << cache_.osdJacobian.block(index * jacobianDim_, 0, jacobianDim_, getDof())
                     * cache_.dcJacobianInvs[it->second.second]
              << std::endl;
  }
  // (2) Check the OSD dynamics:
  std::cout << "The OSD dynamics equation is: " << std::endl
            << cache_.lambdaMatrix * cache_.osdAcc + cache_.rhoOne + cache_.rhoTwo - cache_.osdTau << std::endl;

  std::cout << "The OSD dynamics equation two is: " << std::endl
            << cache_.lambdaMatrix * cache_.osdAcc - cache_.osdTau << std::endl;

  std::cout << "The OSD Acc is: " << std::endl << cache_.osdAcc << std::endl;

  std::cout << "The OSD dynamics equation force is: " << std::endl << cache_.osdTau << std::endl;
}
