#include "mi_impactPredictor.h"

mi_impactPredictor::mi_impactPredictor(mc_rbdyn::Robot & robot,
                                       // std::shared_ptr<rbd::ForwardDynamics> & fdPtr,
                                       std::string impactBodyName,
                                       bool linearJacobian,
                                       double impactDuration,
                                       double coeRes)
: robot_(robot), linearJacobian_(linearJacobian), impactDuration_(impactDuration), coeRes_(coeRes)
{

  std::cout << "The impact predictor constuctor is started." << std::endl;
  setImpactBody(impactBodyName);
  std::cout << "The impact body name is: " << getImpactBody_() << std::endl;
  // osdPtr_ = std::make_shared<mi_osd>(getRobot(), fdPtr, useLinearJacobian_());
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

void mi_impactPredictor::run(const Eigen::Vector3d & surfaceNormal)
{
  if(cache_.grfContainer.size() < 3)
  {
    throw std::runtime_error("Impact-predictor: There are too less end-effector defined");
  }
  assert(cache_.grfContainer.size() == getOsd_()->getEeNum());

  // Update the equations of motions
  osdPtr_->update();
  // std::cout << "OSD updated. " << std::endl;
  Eigen::Matrix3d tempProjector = surfaceNormal * surfaceNormal.transpose();

  Eigen::VectorXd alpha = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alphaD);
  const auto & impactBodyValuesPtr = cache_.grfContainer.find(getImpactBody_());

  //(0.0) update impact body-velocity jump

  impactBodyValuesPtr->second.deltaV = -(getCoeRes_() + 1) * tempProjector * getOsd_()->getJacobian(getImpactBody_())
                                       * (alpha + alphaD * getImpactDuration_());
  //(0.1) update impact body-velocity impulsive force
  impactBodyValuesPtr->second.impulseForce = (1 / getImpactDuration_())
                                             * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_())
                                             * getOsd_()->getDcJacobianInv(getImpactBody_()) * getEeVelocityJump();

  //(0.2) update impact body-velocity induced joint velocity jump
  impactBodyValuesPtr->second.deltaQDot = getOsd_()->getDcJacobianInv(getImpactBody_()) * getEeVelocityJump();

  // add the impact body joint velocity first
  cache_.qVelJump = impactBodyValuesPtr->second.deltaQDot;

  //(0.3) update impact body-velocity induced impulsive joint torque
  // * Update the impulsive force of end-effectors with established contact
  impactBodyValuesPtr->second.deltaTau = getOsd_()->getJacobian(getImpactBody_()).transpose() * getImpulsiveForce();

  // add the impact body impulsive force first
  cache_.tauJump = impactBodyValuesPtr->second.deltaTau;

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

    //(1.0) update impact body-velocity induced end-effector velocity jump and joint velocity jump
    /*
        it->second.deltaV = getOsd_()->getLambdaMatrixInv().block(getOsd_()->nameToIndex_(it->first) * dim,
                                                                  getOsd_()->nameToIndex_(getImpactBody_()) * dim, dim,
       dim)
                            * getImpulsiveForce();
          */
    it->second.deltaV = getOsd_()->getLambdaMatrixInv(it->first, getImpactBody_()) * getImpulsiveForce();

    it->second.deltaQDot = getOsd_()->getDcJacobianInv(it->first) * it->second.deltaV;
    //(1.1) update impact body-velocity induced end-effector impulsive force
    // End-effector reaction force:

    // Add to dq
    cache_.qVelJump += it->second.deltaQDot;
    cache_.tauJump += it->second.deltaTau;
  } // End of looping through the container

  // (2.0) Calculate the propagated impulse
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    if(it->first == getImpactBody_())
    {
      // We skip the impact body
      continue;
    }
    it->second.impulseForce.setZero();
    if(it->second.contact())
    { // If this body has an established contact
      // loop over the entire row according to the OSD model
      for(auto idx = cache_.grfContainer.begin(); idx != cache_.grfContainer.end(); ++idx)
      {
        it->second.impulseForce += getOsd_()->getLambdaMatrix(it->first, idx->first) * getEeVelocityJump(idx->first);
      } // end of row calculation
      it->second.impulseForce = (1 / getImpactDuration_()) * it->second.impulseForce;
      it->second.deltaTau = getOsd_()->getJacobian(it->first).transpose() * it->second.impulseForce;
    } // end of contact force jump

  } // end of impulsive force and torque loop

  // (3.0) update the jacobians for mc_rtc

  Eigen::MatrixXd tempJDeltaAlpha;
  tempJDeltaAlpha.resize(getOsd_()->getDof(), getOsd_()->getJacobianDim());
  tempJDeltaAlpha.setZero();

  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    if(it->first == getImpactBody_())
    {
      // We skip the impact body
      continue;
    }
    tempJDeltaAlpha +=
        getOsd_()->getDcJacobianInv(it->first) * getOsd_()->getLambdaMatrixInv(it->first, getImpactBody_());

  } // end of temp J delta alpha

  Eigen::MatrixXd tempJDeltaTau;
  tempJDeltaTau.resize(getOsd_()->getDof(), 3);
  tempJDeltaTau.setZero();

  //for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  //{
  auto it = cache_.grfContainer.find("l_sole");
 /* 
    if(it->first == getImpactBody_())
    {
      // Skip impact body not in contact
      continue;
    }
  */  
    Eigen::Matrix3d temp;
    temp.setZero();

    for(auto idx = cache_.grfContainer.begin(); idx != cache_.grfContainer.end(); ++idx)
    {
/*	    
      if(it->first == getImpactBody_())
      {
        // Skip body not in contact
        continue;
      }
 */     
      temp += getOsd_()->getLambdaMatrix(it->first, idx->first)
              * getOsd_()->getLambdaMatrixInv(idx->first, getImpactBody_());

    } // end of inner loop
    tempJDeltaTau += getOsd_()->getJacobian(it->first).transpose() * temp;
 // } // end of temp J delta tau
  

  cache_.jacobianDeltaAlpha =
      -(1 + getCoeRes_())
      * (getOsd_()->getDcJacobianInv(getImpactBody_())
         + (1 / getImpactDuration_()) * tempJDeltaAlpha * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_()) * getOsd_()->getDcJacobianInv(getImpactBody_())
	 )
      * tempProjector * getOsd_()->getJacobian(getImpactBody_());
      

  
  cache_.jacobianDeltaTau = 
	  -(1 + getCoeRes_())
                            * (1 / getImpactDuration_()) * tempJDeltaTau
                            * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_()) 
			    * getOsd_()->getDcJacobianInv(getImpactBody_())
			    * tempProjector
                            * getOsd_()->getJacobian(getImpactBody_());
			    
 /* 
  cache_.jacobianDeltaTau = 
	  -(1 + getCoeRes_())
                            * (getOsd_()->getJacobian(getImpactBody_()).transpose() 
			       +(1 / getImpactDuration_()) * tempJDeltaTau
			      )
                            * getOsd_()->getEffectiveLambdaMatrix(getImpactBody_()) 
			    * getOsd_()->getDcJacobianInv(getImpactBody_())
			    * tempProjector
                            * getOsd_()->getJacobian(getImpactBody_());
			    
*/			    
}
