#include "mi_osd.h"

mi_osd::mi_osd(// const dart::dynamics::SkeletonPtr & robotPtr,
		mc_rbdyn::Robot & robot, 
		bool linearJacobian) : robot_(robot)//robotPtr_(robotPtr), 
{
  std::cout << "The osd dynamics constructor is called " << std::endl;
  linearJacobian_ = linearJacobian;
  // Initilize the forward dynamics:
  FDPtr_ = std::make_shared<rbd::ForwardDynamics>(getRobot().mb());

  std::cout << "The FD constructor is built." << std::endl;
  // Create a local copy to avoid touching the mc_rtc controller robot. 
  //rbd::MultiBodyConfig & tempMbc = getRobot().mbc();
  rbd::forwardKinematics(getRobot().mb(), getRobot().mbc() );
  rbd::forwardVelocity(getRobot().mb(),getRobot().mbc() );
  rbd::forwardAcceleration(getRobot().mb(),getRobot().mbc() );
  FDPtr_->forwardDynamics(getRobot().mb(),getRobot().mbc() );
  //FDPtr_->computeH(getRobot().mb(), getRobot().mbc());
  std::cout << "The masss matrix is built." << std::endl;
  // Initialize the Jacobians
  int mRows = static_cast<int>(getFD()->H().rows());
  int mCols = static_cast<int>(getFD()->H().cols());

  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_M(getFD()->H());
  cache_.invMassMatrix.resize(mRows, mCols);
  cache_.invMassMatrix = lu_decomp_M.inverse();

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

  std::cout << "The stack of Jacobians are about to be built." << std::endl;
  // I can use a configuration file to specify the end-effectors later.
  /*
  cache_.jacobians.insert(
      std::make_pair(
        strdup("l_ankle"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "l_ankle" ),
          0
          )
        )
      );

  cache_.jacobians.insert(
      std::make_pair(
        strdup("r_ankle"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "r_ankle" ),
          1
          )
        )
      );

  cache_.jacobians.insert(
      std::make_pair(
        strdup("l_wrist"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "l_wrist" ),
          2
          )
        )
      );
  cache_.jacobians.insert(
      std::make_pair(
        strdup("r_wrist"),
        std::make_pair(
          std::make_shared<rbd::Jacobian>(getRobot()->mb, "r_wrist" ),
          3
          )
        )
      );
*/
  std::string l_ankle_string("l_ankle");

  cache_.jacobians[l_ankle_string] = std::make_pair(std::make_shared<rbd::Jacobian>(getRobot().mb(), "r_wrist"), 0);

  std::string r_ankle_string("r_ankle");
  cache_.jacobians[r_ankle_string] = std::make_pair(std::make_shared<rbd::Jacobian>(getRobot().mb(), "r_wrist"), 1);

  std::string l_wrist_string("l_wrist");
  cache_.jacobians[l_wrist_string] = std::make_pair(std::make_shared<rbd::Jacobian>(getRobot().mb(), "r_wrist"), 2);

  std::string r_wrist_string("r_wrist");
  cache_.jacobians[r_wrist_string] = std::make_pair(std::make_shared<rbd::Jacobian>(getRobot().mb(), "r_wrist"), 3);

  std::cout << "The stack of Jacobians are built." << std::endl;

  eeNum_ = static_cast<int>(cache_.jacobians.size());

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
void mi_osd::updateCache_()
{
  // Read from the robot:
  std::cout << "Updating OSD cache..." << std::endl;

  // Update the mass matrix inverse
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_M(getFD()->H());
  cache_.invMassMatrix = lu_decomp_M.inverse();
  std::cout << "Inverse mass matrix calculated, there are row: " << getInvMassMatrix().rows()
            << ", col: " << getInvMassMatrix().cols() << std::endl;

  for(auto it = cache_.jacobians.begin(); it != cache_.jacobians.end(); ++it)
  {
    int ii = it->second.second;
    // std::cout << it->first->getName() << " has a local index: " << it->second << std::endl;
    // cache_.osdJacobian.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) = it->first->getJacobian();
    //
    auto tempJacobian = it->second.first->jacobian(getRobot().mb(), getRobot().mbc());
    auto tempJacobianDot = it->second.first->jacobianDot(getRobot().mb(), getRobot().mbc());
    // std::cout<<"tempJacobian size is: "<<tempJacobian.rows()<<", "<<tempJacobian.cols()<<std::endl;
    std::cout << "Jacobian matrix calculated..." << std::endl;
    Eigen::MatrixXd tempFullJacobian, tempFullJacobianDot;
    tempFullJacobian.resize(jacobianDim_, getDof());
    tempFullJacobianDot.resize(jacobianDim_, getDof());

    // std::cout<<"temp Full Jacobian size is: "<<tempFullJacobian.rows()<<", "<<tempFullJacobian.cols()<<std::endl;
    if(useLinearJacobian_())
    {
      it->second.first->fullJacobian(getRobot().mb(), 
		      tempJacobian.block(3, 0, 3, tempJacobian.cols()),
                                     tempFullJacobian);
      it->second.first->fullJacobian(getRobot().mb(), 
		      tempJacobian.block(3, 0, 3, tempJacobianDot.cols()),
                                     tempFullJacobianDot);
      /*
      cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) = 
	      getDartRobot()->getBodyNode(it->first)->getCOMLinearAcceleration();
	      */

      cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) = 
	    (getRobot().mbc().bodyPosW[getRobot().mb().bodyIndexByName(it->first)]
	    *getRobot().mbc().bodyAccB
	    [
	    getRobot().mb().bodyIndexByName(it->first)
	    ].vector()).linear();

      cache_.osdVel.segment(ii * jacobianDim_, jacobianDim_) = 
	    getRobot().mbc().bodyVelW
	    [
	    getRobot().mb().bodyIndexByName(it->first)
	    ].linear();

    }
    else
    {
      it->second.first->fullJacobian(getRobot().mb(), tempJacobian, tempFullJacobian);
      it->second.first->fullJacobian(getRobot().mb(), tempJacobianDot, tempFullJacobianDot);
      cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) = 
	    (getRobot().mbc().bodyPosW[getRobot().mb().bodyIndexByName(it->first)]
	    *getRobot().mbc().bodyAccB
	    [
	    getRobot().mb().bodyIndexByName(it->first)
	    ].vector()).vector();

/*
      cache_.osdAcc.segment(ii * jacobianDim_, jacobianDim_) = 
	    getDartRobot()->getBodyNode(it->first)->getCOMSpatialAcceleration();

*/
      
      cache_.osdVel.segment(ii * jacobianDim_, jacobianDim_) = 
			      getRobot().mbc().bodyVelW
			      [
			      getRobot().mb().bodyIndexByName(it->first)
			      ].vector();


    }

    std::cout << "Full Jacobian matrix calculated..." << std::endl;

    std::cout << "" << std::endl;
    cache_.osdJacobian.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) = tempFullJacobian;
    cache_.osdJacobianDot.block(ii * jacobianDim_, 0, jacobianDim_, getDof()) = tempFullJacobianDot;
   /* 
    cache_.osdTau.block(ii * jacobianDim_, 0, jacobianDim_, 1) = 
	    getDartRobot()->getForce(
			    getDartBodyIndex_(it->first)
			    );
			    */

        
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_J_temp(
        cache_.osdJacobian.block(ii * jacobianDim_, 0, jacobianDim_, getDof()));
    std::cout << it->first << " Jacobian has rank: " << lu_decomp_J_temp.rank() << std::endl;

  }

  std::cout << "OSD Jacobian calculated..." << std::endl;

  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_J(cache_.osdJacobian);
  std::cout << "The stacked OSD Jacobian has row: " << cache_.osdJacobian.rows()
            << ", col: " << cache_.osdJacobian.cols() << ", rank: " << lu_decomp_J.rank() << std::endl;

  cache_.lambdaMatrixInv = cache_.osdJacobian * getInvMassMatrix() * (cache_.osdJacobian.transpose());

  std::cout << "Inverse Lambda mass matrix calculated..." << std::endl;
  // Compute the Lambda matrix component-wise
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
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_lambda(cache_.lambdaMatrix);
  std::cout << "The OSD mass matrix has row: " << cache_.lambdaMatrix.rows() << ", col: " << cache_.lambdaMatrix.cols()
            << ", rank: " << lu_decomp_lambda.rank() << std::endl;
  // std::cout<<"The OSD MassMatrix is: "<< std::endl<<cache_.lambdaMatrix<<std::endl;
  auto tempJointTorque = rbd::dofToVector(
			    getRobot().mb(), 
			    getRobot().mbc().jointTorque
			    );
  std::cout<<"The joint torque is: "<< std::endl<<tempJointTorque<<std::endl;
  std::cout<<"The number of ee is: "<<getEeNum()<<std::endl;
  // Update the dynamically consistent Jacobian inverse:
  for(int ii = 0; ii < getEeNum(); ii++)
  {
    //std::cout<<"ii: "<<ii<<std::endl;
    cache_.effectiveLambdaMatrices[ii] =
        cache_.lambdaMatrix.block(ii * jacobianDim_, 0, jacobianDim_, getEeNum() * jacobianDim_) * cache_.osdJacobian;

    cache_.dcJacobianInvs[ii] = (cache_.effectiveLambdaMatrices[ii] * getInvMassMatrix()).transpose();

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_dcJ(cache_.dcJacobianInvs[ii]);
        // Update the OSD components: 
   cache_.rhoTwo.segment(ii * jacobianDim_, jacobianDim_) = 
	   cache_.dcJacobianInvs[ii].transpose()
	   *getFD()->C();

   cache_.osdTau.segment(ii * jacobianDim_, jacobianDim_) = 
	   cache_.dcJacobianInvs[ii].transpose()
  	  *tempJointTorque;
   std::cout << "The dynamically consistent Jacobian " << ii << " has row: " << cache_.dcJacobianInvs[ii].rows()
	   << ", col: " << cache_.dcJacobianInvs[ii].cols() << ", rank: " << lu_decomp_dcJ.rank() << std::endl;

  }
 // Update the rest of the OSD components 

  std::cout<<"Rho One is to be calculated.";

  
  Eigen::VectorXd jointVel = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  cache_.rhoOne = - cache_.lambdaMatrix
	  *cache_.osdJacobianDot
	  *jointVel;

  std::cout<<"Rho One is calculated.";
  
/*
  cache_.rhoTwo = cache_.dcJacobianInvs.transpose()
  	  *getFD()->C();

  cache_.osdF = cache_.dcJacobianInvs.transpose()
  	  *cache_.osdTau;

*/



 // (0) This is a temporary test: 

  std::cout<<"Inverse of the OSD inertia matrix is: "<<std::endl<<cache_.lambdaMatrixInv<<std::endl;
  std::cout<<"The OSD inertia matrix is: "<<std::endl<<cache_.lambdaMatrix<<std::endl;
  std::cout<<"Multipulication of Lambda*Lambda_inv is: "<<std::endl<<cache_.lambdaMatrix*cache_.lambdaMatrixInv<<std::endl;
  // (1) dc Jacobian inverse multiplies J should be close to identity? 
  
  for(auto it = cache_.jacobians.begin(); it != cache_.jacobians.end(); ++it){
  
    int index = it->second.second;
    // J*dc_J_inv:
    std::cout<<"Test "<<it->first<<" Jacobian multiplies DC Jacobian inverse: "<<std::endl<<cache_.osdJacobian.block(index * jacobianDim_, 0, jacobianDim_, getDof())
	    *cache_.dcJacobianInvs[it->second.second]
	    <<std::endl;
  }
  // (2) Check the OSD dynamics: 
  std::cout<<"The OSD dynamics equation is: "<<std::endl<<
  cache_.lambdaMatrix*cache_.osdAcc + cache_.rhoOne + cache_.rhoTwo - cache_.osdTau<<std::endl;

  std::cout<<"The OSD dynamics equation two is: "<<std::endl<<
  cache_.lambdaMatrix*cache_.osdAcc - cache_.osdTau<<std::endl;

  std::cout<<"The OSD Acc is: "<<std::endl<<
   cache_.osdAcc<<std::endl;
  std::cout<<"The OSD Vel is: "<<std::endl<<
   cache_.osdVel<<std::endl;



  std::cout<<"The OSD dynamics equation force is: "<<std::endl<<
   cache_.osdTau<<std::endl;



}
