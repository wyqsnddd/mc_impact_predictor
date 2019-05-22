#include "mi_impactPredictor.h"

mi_impactPredictor::mi_impactPredictor(mc_rbdyn::Robot & robot,
                                       const std::string & impactBodyName,
                                       bool linearJacobian,
                                       double impactDuration,
                                       double coeRes)
: robot_(robot), linearJacobian_(linearJacobian), impactDuration_(impactDuration), coeRes_(coeRes)
{

  std::cout << "The impact predictor constuctor is started." << std::endl;
  setImpactBody(impactBodyName);
  // osdPtr_ = std::make_shared<mi_osd>(robotPtr, getRobot(), useLinearJacobian_());
  osdPtr_ = std::make_shared<mi_osd>(getRobot(), useLinearJacobian_());
  std::cout << "The impact predictor constuctor is finished." << std::endl;
  std::cout << "The impact duration is: " << getImpactDuration_() << ", the coeres is: " << getCoeRes_() << std::endl;

  cache_.newLeeImpulse.resize(getOsd_()->getJacobianDim());
  cache_.newLeeImpulse.setZero();
  cache_.newReeImpulse.resize(getOsd_()->getJacobianDim());
  cache_.newReeImpulse.setZero();
  cache_.new_eeLeeImpulse.resize(getOsd_()->getJacobianDim());
  cache_.new_eeLeeImpulse.setZero();
  cache_.new_eeReeImpulse.resize(getOsd_()->getJacobianDim());
  cache_.new_eeReeImpulse.setZero();

  cache_.eeImpulse.resize(getOsd_()->getJacobianDim());
  cache_.eeImpulse.setZero();
  cache_.eeVelJump.resize(getOsd_()->getJacobianDim());
  cache_.eeVelJump.setZero();

  cache_.qVelJump.resize(getOsd_()->getDof());
  cache_.qVelJump.setZero();
  cache_.tauJump.resize(getOsd_()->getDof());
  cache_.tauJump.setZero();
}
void mi_impactPredictor::resetDataStructure()
{
  getOsd_()->resetDataStructure();
  cache_.eeVelJump.setZero();
  cache_.qVelJump.setZero();
  cache_.tauJump.setZero();
  cache_.eeImpulse.setZero();
  cache_.newLeeImpulse.setZero();
  cache_.newReeImpulse.setZero();
  cache_.new_eeLeeImpulse.setZero();
  cache_.new_eeReeImpulse.setZero();

  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it){
	  it->second.deltaV.setZero();
	  it->second.impulseForce.setZero();
	  it->second.accForce.setZero();
  }
}

void mi_impactPredictor::initializeDataStructure(int numEE)
{
  getOsd_()->initializeDataStructure(numEE);
}

void mi_impactPredictor::run()
{
  if(cache_.grfContainer.size() < 3)
  {
    throw std::runtime_error("Impact-predictor: There are too less end-effector defined");
  }
  assert(cache_.grfContainer.size() == getOsd_->getEeNum());

  // Update the equations of motions
  std::cout << "mi_impactPredictor::update() is called. " << std::endl;
  osdPtr_->update();
  // Update the data in the cache
  // * Update the end-effector velocityJump
  std::cout << "OSD updated. " << std::endl;
  Eigen::VectorXd alpha = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alphaD);
  std::cout << "q_d is: " << alpha.transpose() << std::endl;
  std::cout << "The velocity of " << getImpactBody_() << " is: " << getRobot().bodyVelW(getImpactBody_()) << std::endl;
  // std::cout<<"q_dd is: "<<alphaD.transpose()<<std::endl;
  // std::cout<<" The Jacobian is: "<<getOsd_()->getJacobian(getImpactBody_())<<std::endl;

  cache_.eeVelJump =
      -(getCoeRes_() + 1) * getOsd_()->getJacobian(getImpactBody_()) * (alpha + alphaD * getImpactDuration_());
/*
  cache_.eeVelJump =
      -(getCoeRes_() + 1) * getOsd_()->getJacobian(getImpactBody_()).block(1, 0, 1, getOsd_()->getJacobianDim())
      * (alpha + alphaD * getImpactDuration_());
*/

  cache_.eeImpulse = (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_())
                     * getOsd_()->getDcJacobianInv(getImpactBody_()) * cache_.eeVelJump;

  // std::cout<<"The impact body is: "<<getImpactBody_()<<std::endl;
  // * Update the joint velocity jump
  cache_.qVelJump = getOsd_()->getDcJacobianInv(getImpactBody_()) * cache_.eeVelJump;
  // cache_.qVelJump = getOsd_()->getNewDcJacobianInv(getImpactBody_()) * cache_.eeVelJump;
  std::cout << "The predicted impact body velocity jump is: " << std::endl << cache_.eeVelJump.transpose() << std::endl;
  // * Update the impulsive force

  /*
  cache_.eeImpulse = (1 / getImpactDuration_())
    *getOsd_()->getEffectiveLambdaMatrix(getImpactBody_())
    * getOsd_()->getInvMassMatrix().transpose()
    * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_()).transpose()
    * cache_.eeVelJump;
    */
  /*
   cache_.eeImpulse = (1 / getImpactDuration_())
     * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_())
     * getOsd_()->getDcJacobianInv(getImpactBody_())
     * cache_.eeVelJump;
 */

  cache_.tauJump = getOsd_()->getJacobian(getImpactBody_()).transpose() * cache_.eeImpulse;
  // * Update the impulsive force of end-effectors with established contact

  // std::cout<<"The predicted joint velocity jump is: "<<std::endl<<cache_.qVelJump.transpose()<<std::endl;
  std::cout << "The predicted end-effector impulsive force is: " << std::endl
            << cache_.eeImpulse.transpose() << std::endl;

  std::cout << "------------------Impact body impulsive forces ------------------------------------" << std::endl;
  // Update the ground reaction forces:
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {

    // End-effector velocity jump:
    it->second.deltaV = getOsd_()->getJacobian(it->first) * cache_.qVelJump;
    // End-effector reaction force:
    it->second.impulseForce = (1 / getImpactDuration_())
                              * getOsd_()->getEffectiveLambdaMatrix(it->first)
                              // cache_.eeVelJump;
                              * getOsd_()->getDcJacobianInv(getImpactBody_()) * cache_.eeVelJump;

    Eigen::Vector3d tempGRF_two =
        (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix(it->first) * cache_.qVelJump;

    Eigen::Vector3d tempGRF_three = getOsd_()->getDcJacobianInv(it->first).transpose() * cache_.tauJump;

    std::cout << "The predicted GRF impulsive force of " << it->first << " is: " << it->second.impulseForce.transpose()
              << std::endl
              << " velocity jump is: " << it->second.deltaV.transpose() << std::endl
              << "The predicted GRF impulsive force two is " << tempGRF_two.transpose() << std::endl
              << "The predicted GRF impulsive force three is " << tempGRF_three.transpose() << std::endl;
  }
  // tempTest_();
  tempTestBody_();
  tempTestEe_();
  tempTestAcc_();
  // tempTestAccEe_();
}
void mi_impactPredictor::tempTestAccEe_()
{
  Eigen::VectorXd q_dot = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  Eigen::VectorXd q_ddot = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alphaD);
  int dim = getOsd_()->getJacobianDim();

  Eigen::VectorXd osdAcc =
      getOsd_()->getJacobian(getImpactBody_()) * q_ddot + getOsd_()->getJacobianDot(getImpactBody_()) * q_dot;
  /*
    Eigen::VectorXd osdForce = getOsd_()->getLambdaMatrix().block(
        0,
       getOsd_()->nameToIndex_(getImpactBody_())*dim,
       dim*getOsd_()->getEeNum(),
       dim
        )
      *osdAcc;
  */
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    /*
  it->second.accForce = osdForce.segment(
      getOsd_()->nameToIndex_(it->first)*dim,
      dim
      );
  }
  */

    it->second.accForce =
        getOsd_()->getLambdaMatrix(getOsd_()->nameToIndex_(it->first), getOsd_()->nameToIndex_(getImpactBody_()))
        * osdAcc;
  }
}

void mi_impactPredictor::tempTestAcc_()
{
  Eigen::VectorXd q_dot = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  Eigen::VectorXd q_ddot = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alphaD);
  int dim = getOsd_()->getJacobianDim();
  Eigen::VectorXd osdAcc;
  osdAcc.resize(dim * getOsd_()->getEeNum());
  osdAcc.setZero();
  
    osdAcc.segment(
         getOsd_()->nameToIndex_("l_sole")*dim,
         dim
       ) =
       getOsd_()->getJacobian("l_sole")*q_ddot +
       getOsd_()->getJacobianDot("l_sole")*q_dot;
  
/*
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    // Calculate the acceleration for all the end effectors
    osdAcc.segment(getOsd_()->nameToIndex_(it->first) * dim, dim) =
        getOsd_()->getJacobian(it->first) * q_ddot + getOsd_()->getJacobianDot(it->first) * q_dot;
  }
*/
  Eigen::VectorXd osdForce = getOsd_()->getLambdaMatrix() * osdAcc;
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {

    it->second.accForce = osdForce.segment(getOsd_()->nameToIndex_(it->first) * dim, dim);
  }
}
void mi_impactPredictor::tempTestEe_()
{

  // For the established contacts:
  int dim = getOsd_()->getJacobianDim();

  Eigen::VectorXd delta_vel_l_ankle =
      getOsd_()->getLambdaMatrixInv().block(getOsd_()->nameToIndex_("l_sole") * dim,
                                            getOsd_()->nameToIndex_(getImpactBody_()) * dim, dim, dim)
      * cache_.eeImpulse;

  Eigen::VectorXd delta_vel_r_ankle =
      getOsd_()->getLambdaMatrixInv().block(getOsd_()->nameToIndex_("r_sole") * dim,
                                            getOsd_()->nameToIndex_(getImpactBody_()) * dim, dim, dim)
      * cache_.eeImpulse;

  cache_.new_eeLeeImpulse = (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix("l_sole")
                            * getOsd_()->getDcJacobianInv("l_sole") * delta_vel_l_ankle;

  cache_.new_eeReeImpulse = (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix("r_sole")
                            * getOsd_()->getDcJacobianInv("r_sole") * delta_vel_r_ankle;
}

void mi_impactPredictor::tempTestBody_()
{

  Eigen::VectorXd bodyImpulseForce = (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix("body")
                                     * getOsd_()->getDcJacobianInv(getImpactBody_()) * cache_.eeVelJump;

  // For the established contacts:
  int dim = getOsd_()->getJacobianDim();

  Eigen::VectorXd delta_vel_l_ankle =
      getOsd_()->getLambdaMatrixInv().block(getOsd_()->nameToIndex_("l_sole") * dim,
                                            getOsd_()->nameToIndex_("body") * dim, dim, dim)
      * bodyImpulseForce;

  Eigen::VectorXd delta_vel_r_ankle =
      getOsd_()->getLambdaMatrixInv().block(getOsd_()->nameToIndex_("r_sole") * dim,
                                            getOsd_()->nameToIndex_("body") * dim, dim, dim)
      * bodyImpulseForce;

  cache_.newLeeImpulse = (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix("l_sole")
                         * getOsd_()->getDcJacobianInv("l_sole") * delta_vel_l_ankle;

  cache_.newReeImpulse = (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix("r_sole")
                         * getOsd_()->getDcJacobianInv("r_sole") * delta_vel_r_ankle;
}
void mi_impactPredictor::tempTest_()
{
  /*
  Eigen::VectorXd tempEeAcc;
  tempEeAcc.resize(getOsd_()->getEeNum()*getOsd_()->getJacobianDim());
/// This might be a bug in the future, basically I assumed that I know the order
  for (int ii=0; ii<getOsd_()->getEeNum(); ii++){
    tempEeAcc.segment(ii*getOsd_()->getJacobianDim(), getOsd_()->getJacobianDim()) =

}
*/
  // Test the end-effector induced ground reaction forces.
  auto tempImpactBodyAcceleration = getRobot().mbc().bodyPosW[getRobot().mb().bodyIndexByName(getImpactBody_())]
                                    * getRobot().mbc().bodyAccB[getRobot().mb().bodyIndexByName(getImpactBody_())];
  std::cout << "The impact body acceleration is: " << tempImpactBodyAcceleration.linear().transpose() << std::endl;
  std::cout << "------------------Impact body Acc ------------------------------------" << std::endl;
  // Note that we need to deduct the gravity force.
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    std::cout << "The pseudo mass matrix between " << it->first << " and " << getImpactBody_() << " is: " << std::endl
              << getOsd_()->getLambdaMatrix(getImpactBody_(), it->first) << std::endl;

    Eigen::Vector3d tempGRF = getOsd_()->getDcJacobianInv(it->first).transpose() * getOsd_()->getFD()->H()
                              * getOsd_()->getDcJacobianInv(getImpactBody_())
                              // getOsd_()->getLambdaMatrix(getImpactBody_(), it->first)
                              // getOsd_()->getLambdaMatrix(it->first, getImpactBody_())
                              // getOsd_()->getCrossLambdaMatrix( it->first, getImpactBody_())
                              // getOsd_()->getCrossLambdaMatrix( getImpactBody_(), it->first)
                              * tempImpactBodyAcceleration.linear();

    /*
     // (J_dc_inv_i_t M J_dc_inv_m)
     Eigen:: Vector3d tempGRF =
         (getOsd_()->getDcJacobianInv(it->first).transpose()
     * getOsd_()->getFD()->H()//getInvMassMatrix()
     *getOsd_()->getDcJacobianInv(getImpactBody_()))
     * tempImpactBodyAcceleration.linear();
   */
    /*
     // (J_m M_inv Ji_transpose).inv
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_lambda_component(
            getOsd_()->getJacobian(getImpactBody_())
            *getOsd_()->getInvMassMatrix()
            *getOsd_()->getJacobian(it->first).transpose()
            );
     Eigen:: Vector3d tempGRF =
       lu_decomp_lambda_component.inverse()
     * tempImpactBodyAcceleration.linear();
 */
    /*
 Eigen::VectorXd tempAccForce = getOsd_()->getLambdaMatrix().block(
     getOsd_()->nameToIndex_(it->first)*getOsd_()->getJacobianDim(),
     0,
     getOsd_()->getJacobianDim(),
     getOsd_()->getJacobianDim()*getOsd_()->getEeNum())
   *getOsd_()->getOsdAcc();
   */
    // it->second.accForce = tempAccForce;

    it->second.accForce = tempGRF;
    /*
        std::cout << "The predicted GRF force of body " << it->first << " due to impact-body Acc of: " <<
       getImpactBody_()
                  << " is: " << std::endl
                  << tempGRF.transpose() << std::endl;
            */
  }
  /*
    std::cout << "------------------Impact body Inverse Acc ------------------------------------" << std::endl;
    // robot.forceSensor("RightFootForceSensor");
    Eigen::VectorXd exForce;
    exForce.resize(getOsd_()->getJacobianDim() * getOsd_()->getEeNum());
    exForce.setZero();

    int lAnkleIndex = getOsd_()->nameToIndex_("l_ankle");
    exForce.segment(lAnkleIndex * getOsd_()->getJacobianDim(), getOsd_()->getJacobianDim()) =
        getRobot().forceSensor("LeftFootForceSensor").force();

    int rAnkleIndex = getOsd_()->nameToIndex_("r_ankle");
    exForce.segment(rAnkleIndex * getOsd_()->getJacobianDim(), getOsd_()->getJacobianDim()) =
        getRobot().forceSensor("RightFootForceSensor").force();

    std::cout << "The real ex force is: " << std::endl << exForce.transpose() << std::endl;
    std::cout << "The predicted ex force is: " << std::endl
              << (getOsd_()->getLambdaMatrix() * getOsd_()->getOsdAcc()).transpose() << std::endl;

    Eigen::VectorXd predictedAcc = getOsd_()->getLambdaMatrixInv() * exForce;

    std::cout << "The predicted Acc is: " << std::endl << predictedAcc.transpose() << std::endl;
    std::cout << "The actual Acc is: " << std::endl << getOsd_()->getOsdAcc().transpose() << std::endl;
    */
}
