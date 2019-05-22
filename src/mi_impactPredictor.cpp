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
  cache_.ini(getOsd_()->getJacobianDim());
}
void mi_impactPredictor::resetDataStructure()
{
  getOsd_()->resetDataStructure();
  cache_.reset();
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
  tempTestEe_();
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

