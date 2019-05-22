#include "mi_impactPredictor.h"

mi_impactPredictor::mi_impactPredictor(mc_rbdyn::Robot & robot,
                                       std::string impactBodyName,
                                       bool linearJacobian,
                                       double impactDuration,
                                       double coeRes)
: robot_(robot), linearJacobian_(linearJacobian), impactDuration_(impactDuration), coeRes_(coeRes)
{

  std::cout << "The impact predictor constuctor is started." << std::endl;
  setImpactBody(impactBodyName);
  std::cout << "The impact body name is: " << getImpactBody_() << std::endl;
  osdPtr_ = std::make_shared<mi_osd>(getRobot(), useLinearJacobian_());

  std::cout << "The impact predictor constuctor is finished." << std::endl;
  std::cout << "The impact duration is: " << getImpactDuration_() << ", the coeres is: " << getCoeRes_() << std::endl;
}

bool mi_impactPredictor::addEndeffector(std::string eeName)
{

  unsigned eeNum = static_cast<unsigned>(cache_.grfContainer.size());

  if(useLinearJacobian_())
  {
    cache_.grfContainer[eeName] = {false, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof()),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof())};
  }
  else
  {
    cache_.grfContainer[eeName] = {false, Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof()),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof())};
  }

  if(!getOsd_()->addEndeffector(eeName))
  {
    throw std::runtime_error(std::string("OSd failed to add endeffector!") + eeName);
  }

  unsigned eeNumNew = static_cast<unsigned>(cache_.grfContainer.size());
  if((eeNum == eeNumNew - 1) && (eeNumNew == static_cast<unsigned>(getOsd_()->getEeNum())))
  {
    return true;
  }
  else
  {
    std::cout << "The ee container start with " << eeNum << ", now it has " << eeNumNew << " endeffectors. The osd has "
              << getOsd_()->getEeNum() << " end-effectors. " << std::endl;
    return false;
  }
}

void mi_impactPredictor::run()
{
  if(cache_.grfContainer.size() < 3)
  {
    throw std::runtime_error("Impact-predictor: There are too less end-effector defined");
  }
  assert(cache_.grfContainer.size() == getOsd_()->getEeNum());

  // Update the equations of motions
  osdPtr_->update();
  // std::cout << "OSD updated. " << std::endl;

  Eigen::VectorXd alpha = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alphaD);
  const auto & impactBodyValuesPtr = cache_.grfContainer.find(getImpactBody_());

  impactBodyValuesPtr->second.deltaV =
      -(getCoeRes_() + 1) * getOsd_()->getJacobian(getImpactBody_()) * (alpha + alphaD * getImpactDuration_());

  impactBodyValuesPtr->second.impulseForce = (1 / getImpactDuration_())
                                             * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_())
                                             * getOsd_()->getDcJacobianInv(getImpactBody_()) * getEeVelocityJump();

  impactBodyValuesPtr->second.deltaQDot = getOsd_()->getDcJacobianInv(getImpactBody_()) * getEeVelocityJump();

  // add the impact body joint velocity first
  cache_.qVelJump = impactBodyValuesPtr->second.deltaQDot;
  // cache_.qVelJump = getOsd_()->getNewDcJacobianInv(getImpactBody_()) * cache_.eeVelJump;
  std::cout << "The predicted impact body velocity jump is: " << std::endl << getEeVelocityJump() << std::endl;

  // cache_.tauJump =
  //	  getOsd_()->getJacobian(getImpactBody_()).transpose()
  //	  * cache_.eeImpulse;

  // * Update the impulsive force of end-effectors with established contact
  impactBodyValuesPtr->second.deltaTau = getOsd_()->getJacobian(getImpactBody_()).transpose() * getImpulsiveForce();

  // add the impact body impulsive force first
  cache_.tauJump = impactBodyValuesPtr->second.deltaTau;
  // std::cout<<"The predicted joint velocity jump is: "<<std::endl<<cache_.qVelJump.transpose()<<std::endl;
  std::cout << "The predicted end-effector impulsive force is: " << std::endl
            << getImpulsiveForce() << std::endl; // cache_.eeImpulse.transpose() << std::endl;

  // std::cout << "------------------Impact body impulsive forces ------------------------------------" << std::endl;
  // Update the ground reaction forces:
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    if(it->first == getImpactBody_())
    {
      // We skip the impact body
      continue;
    }

    int dim = getOsd_()->getJacobianDim();
    // End-effector velocity jump:

    it->second.deltaV = getOsd_()->getLambdaMatrixInv().block(getOsd_()->nameToIndex_(it->first) * dim,
                                                              getOsd_()->nameToIndex_(getImpactBody_()) * dim, dim, dim)
                        * getImpulsiveForce();

    // it->second.deltaV = getOsd_()->getJacobian(it->first) * cache_.qVelJump;
    // it->second.deltaV = getOsd_()->getJacobian(it->first) * cache_.qVelJump;
    // End-effector reaction force:
    if(it->second.contact())
    {
      it->second.impulseForce = (1 / getImpactDuration_()) * getOsd_()->getEffectiveLambdaMatrix(it->first)
                                * getOsd_()->getDcJacobianInv(it->first) * getEeVelocityJump(it->first);

      it->second.deltaTau = getOsd_()->getJacobian(it->first).transpose() * it->second.impulseForce;
    }

    // Eigen::Vector3d tempGRF_three = getOsd_()->getDcJacobianInv(it->first).transpose() * cache_.tauJump;

    std::cout << "The predicted GRF impulsive force of " << it->first << " is: " << it->second.impulseForce.transpose()
              << std::endl
              << " velocity jump is: " << it->second.deltaV.transpose() << std::endl;
    //         << "The predicted GRF impulsive force three is " << tempGRF_three.transpose() << std::endl;
    // Add to dq
    cache_.qVelJump += it->second.deltaQDot;
    cache_.tauJump += it->second.deltaTau;

    // Add to dtau
  } // End of looping through the container
  // tempTestEe_();
}
/*
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
*/
