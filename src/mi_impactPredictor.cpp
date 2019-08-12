#include "mi_impactPredictor.h"

mi_impactPredictor::mi_impactPredictor(mc_rbdyn::Robot & robot,
                                       const std::shared_ptr<mi_osd> & osdPtr,
                                       std::string impactBodyName,
                                       bool linearJacobian,
                                       double impactDuration,
                                       double timeStep,
                                       double coeFrictionDeduction,
                                       double coeRes)
: robot_(robot), osdPtr_(osdPtr), linearJacobian_(linearJacobian), impactDuration_(impactDuration), timeStep_(timeStep),
  coeFrictionDeduction_(coeFrictionDeduction), coeRes_(coeRes)
{

  std::cout << "The impact predictor constuctor is started." << std::endl;
  setImpactBody(impactBodyName);
  std::cout << "The impact body name is: " << getImpactBody() << std::endl;
  // osdPtr_ = std::make_shared<mi_osd>(getRobot(), fdPtr, useLinearJacobian_());
  // osdPtr_ = std::make_shared<mi_osd>(getRobot(), useLinearJacobian_());

  std::cout << "The impact predictor constuctor is finished." << std::endl;
  std::cout << "The impact duration is: " << getImpactDuration() << ", the coeres is: " << getCoeRes_() << std::endl;
}
void mi_impactPredictor::setContact(std::string contactBodyName)
{
  const auto & ee = cache_.grfContainer.find(contactBodyName);
  if(ee != (cache_.grfContainer.end()))
  {
    ee->second.setContact();
    cache_.contactEndeffectors.push_back(contactBodyName);
    std::cout << "setContact: " << contactBodyName << ": "<<ee->second.contact() << std::endl;
  }
  else
  {
    std::cout << "setContact: " << contactBodyName << std::endl;
    throw std::runtime_error(std::string("setContact: '") + contactBodyName
                             + std::string("' is not in the container."));
  }
}
const Eigen::VectorXd & mi_impactPredictor::getImpulsiveForce(const std::string & eeName)
{
  const auto & ee = cache_.grfContainer.find(eeName);
  if(ee->second.contact() || (ee->first == getImpactBody()))
  {
    return ee->second.impulseForce;
  }
  else
  {
    throw std::runtime_error(std::string("Predictor: '-") + eeName + std::string("- ' is not in contact."));
  }
}
const sva::ForceVecd & mi_impactPredictor::getImpulsiveForceCOM(const std::string & eeName)
{
  const auto & ee = cache_.grfContainer.find(eeName);
  if(ee->second.contact() || (ee->first == getImpactBody()))
  {
    return ee->second.impulseForceCOM;
  }
  else
  {
    throw std::runtime_error(std::string("Predictor: '-") + eeName + std::string("- ' is not in contact."));
  }
}

const sva::ForceVecd & mi_impactPredictor::getImpulsiveForceCOP(const std::string & eeName)
{
  const auto & ee = cache_.grfContainer.find(eeName);
  if(ee->second.contact())
  {
    return ee->second.impulseForceCOP;
  }
  else
  {
    throw std::runtime_error(std::string("Predictor: '-") + eeName + std::string("- ' is not in contact."));
  }
}
bool mi_impactPredictor::addEndeffector(std::string eeName)
{

  unsigned eeNum = static_cast<unsigned>(cache_.grfContainer.size());

  if(useLinearJacobian_())
  {
    cache_.grfContainer[eeName] = {false,
                                   Eigen::Vector3d::Zero(),
                                   Eigen::Vector3d::Zero(),
                                   sva::ForceVecd(Eigen::Vector6d::Zero()),
                                   sva::ForceVecd(Eigen::Vector6d::Zero()),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof()),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof()),
                                   Eigen::MatrixXd::Zero(3, getOsd_()->getDof())};
  }
  else
  {
    cache_.grfContainer[eeName] = {false,
                                   Eigen::Vector6d::Zero(),
                                   Eigen::Vector6d::Zero(),
                                   sva::ForceVecd(Eigen::Vector6d::Zero()),
                                   sva::ForceVecd(Eigen::Vector6d::Zero()),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof()),
                                   Eigen::VectorXd::Zero(getOsd_()->getDof()),
                                   Eigen::MatrixXd::Zero(6, getOsd_()->getDof())};
  }

  if(!getOsd_()->addEndeffector(eeName))
  {
    throw std::runtime_error(std::string("OSD failed to add endeffector! ") + eeName);
  }

  unsigned eeNumNew = static_cast<unsigned>(cache_.grfContainer.size());
  if(eeNum == eeNumNew - 1)
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
  // osdPtr_->update();
  // std::cout << "OSD updated. " << std::endl;
  Eigen::Matrix3d tempProjector = surfaceNormal * surfaceNormal.transpose();
  Eigen::Matrix3d tempNullProjector = Eigen::Matrix3d::Identity() - tempProjector ;

  //Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes_()) * tempProjector + getCoeFricDe_() * tempNullProjector);
  Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes_()) * tempProjector );
  Eigen::VectorXd alpha = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alphaD);
  const auto & impactBodyValuesPtr = cache_.grfContainer.find(getImpactBody());

  //(0.0) update impact body-velocity jump

  // impactBodyValuesPtr->second.deltaV = -(getCoeRes_() + 1) * tempProjector * getOsd_()->getJacobian(getImpactBody())

  impactBodyValuesPtr->second.deltaV =
      tempReductionProjector * getOsd_()->getJacobian(getImpactBody()) * (alpha + alphaD *getTimeStep());
  //(0.1) update impact body-velocity impulsive force
  impactBodyValuesPtr->second.impulseForce = (1 / getImpactDuration())
                                             * getOsd_()->getEffectiveLambdaMatrix(getImpactBody())
                                             * getOsd_()->getDcJacobianInv(getImpactBody()) * getEeVelocityJump();

  // calculate the equivalent wrench at the COM
  sva::PTransformd X_ee_CoM = sva::PTransformd(getRobot().com()) * getRobot().bodyPosW(getImpactBody()).inv();
  impactBodyValuesPtr->second.impulseForceCOM =
      X_ee_CoM.dualMul(sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce(getImpactBody())));

  //(0.2) update impact body-velocity induced joint velocity jump
  impactBodyValuesPtr->second.deltaQDot = getOsd_()->getDcJacobianInv(getImpactBody()) * getEeVelocityJump();

  // add the impact body joint velocity first
  cache_.qVelJump = impactBodyValuesPtr->second.deltaQDot;

  //(0.3) update impact body-velocity induced impulsive joint torque
  // * Update the impulsive force of end-effectors with established contact
  impactBodyValuesPtr->second.deltaTau = getOsd_()->getJacobian(getImpactBody()).transpose() * getImpulsiveForce();

  // (0.4) Update the Jacobian
  impactBodyValuesPtr->second.jacobianDeltaF = getOsd_()->getEffectiveLambdaMatrix(getImpactBody())
                                               * getOsd_()->getDcJacobianInv(getImpactBody()) * tempReductionProjector
                                               * getOsd_()->getJacobian(getImpactBody());

  // add the impact body impulsive force first
  cache_.tauJump = impactBodyValuesPtr->second.deltaTau;

  // Update the acc force
  /*
  Eigen::Vector3d tempImpactBodyAcceleration =
             getOsd_()->getJacobian(getImpactBody())*alphaD
       + getOsd_()->getJacobianDot(getImpactBody())*alpha;


  impactBodyValuesPtr->second.accForce =
    getOsd_()->getLambdaMatrix(getImpactBody(), getImpactBody())
    *tempImpactBodyAcceleration;
*/

  // std::cout << "------------------Impact body impulsive forces ------------------------------------" << std::endl;
  // Update the ground reaction forces:
  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    if(it->first == getImpactBody())
    {
      // We skip the impact body
      continue;
    }

    //int dim = getOsd_()->getJacobianDim();
    // End-effector velocity jump:

    //(1.0) update impact body-velocity induced end-effector velocity jump and joint velocity jump
    /*
        it->second.deltaV = getOsd_()->getLambdaMatrixInv().block(getOsd_()->nameToIndex_(it->first) * dim,
                                                                  getOsd_()->nameToIndex_(getImpactBody()) * dim, dim,
       dim)
                            * getImpulsiveForce();
          */
    it->second.deltaV = getOsd_()->getLambdaMatrixInv(it->first, getImpactBody()) * getImpulsiveForce();

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
    if(it->first == getImpactBody())
    {
      // We skip the impact body
      continue;
    }
    it->second.impulseForce.setZero();
    // it->second.accForce.setZero();

    // For the body with established contact:
    if(it->second.contact())
    { // If this body has an established contact
      // loop over the entire row according to the OSD model
      for(auto idx = cache_.grfContainer.begin(); idx != cache_.grfContainer.end(); ++idx)
      {
        it->second.impulseForce += getOsd_()->getLambdaMatrix(it->first, idx->first) * getEeVelocityJump(idx->first);
        // it->second.impulseForce += getOsd_()->getLambdaMatrix(it->first, idx->first) *
        // getOsd_()->getLambdaMatrixInv(it->first, getImpactBody())*getImpulsiveForce();

        // std::cout<<"temp is: "<<getOsd_()->getLambdaMatrix( it->first, idx->first) *
        // getOsd_()->getLambdaMatrixInv(idx->first, getImpactBody())<<std::endl;
        /*
          it->second.accForce +=
            getOsd_()->getLambdaMatrix(it->first, idx->first)
            *getRobot().mbc().bodyAccB
            [
            getRobot().mb().bodyIndexByName(idx->first)
            ].linear();
        */
      } // end of row calculation
      it->second.impulseForce = (1 / getImpactDuration()) * it->second.impulseForce;
      it->second.deltaTau = getOsd_()->getJacobian(it->first).transpose() * it->second.impulseForce;

      // calculate the equivalent wrench at the COM
      sva::PTransformd X_ee_CoM = sva::PTransformd(getRobot().com()) * getRobot().bodyPosW(it->first).inv();
      it->second.impulseForceCOM =
          X_ee_CoM.dualMul(sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce(it->first)));

    } // end of contact force jump

    // 2.2 Update the acc force
    /*
     Eigen::Vector3d tempImpactBodyAcceleration =
              getOsd_()->getJacobian(getImpactBody())*alphaD
        + getOsd_()->getJacobianDot(getImpactBody())*alpha;
        */
    /*
          auto tempImpactBodyAcceleration =
                  getRobot().mbc().bodyAccB
                  [
            getRobot().mb().bodyIndexByName(it->first)
                  ];
            */

  } // end of impulsive force and torque loop

  // 2.2 Calculate the sum of the equivalent impulsive forces:
  for(auto ii = cache_.contactEndeffectors.begin(); ii != cache_.contactEndeffectors.end(); ++ii)
  {
    const auto & ee = cache_.grfContainer.find(*ii);

    sva::PTransformd X_0_ee = getRobot().bodyPosW(getImpactBody());
    sva::PTransformd X_0_self = getRobot().bodyPosW(*ii);

    sva::PTransformd X_ee_self = X_0_self * X_0_ee.inv();
    sva::ForceVecd f_ee = X_ee_self.dualMul(sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce()));

    ee->second.impulseForceCOP = sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce(*ii)) + f_ee;

    for(auto jj = cache_.contactEndeffectors.begin(); jj != cache_.contactEndeffectors.end(); ++jj)
    {
      if(jj == ii)
      {
        continue;
      }
      sva::PTransformd X_0_other = getRobot().bodyPosW(*jj);
      sva::PTransformd X_other_self = X_0_self * X_0_other.inv();
      ee->second.impulseForceCOP +=
          X_other_self.dualMul(sva::ForceVecd(Eigen::Vector3d::Zero(), getImpulsiveForce(*jj)));
    }
  }
  // (3.0) update the jacobians for mc_rtc

  // ----------------------------------- Jacobian Alpha
  Eigen::MatrixXd tempJDeltaAlpha;
  tempJDeltaAlpha.resize(getOsd_()->getDof(), getOsd_()->getJacobianDim());
  tempJDeltaAlpha.setZero();

  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    if(it->first == getImpactBody())
    {
      // We skip the impact body
      continue;
    }
    tempJDeltaAlpha +=
        getOsd_()->getDcJacobianInv(it->first) * getOsd_()->getLambdaMatrixInv(it->first, getImpactBody());

  } // end of temp J delta alpha

  cache_.jacobianDeltaAlpha =
      (getOsd_()->getDcJacobianInv(getImpactBody())
       + (1 / getImpactDuration()) * tempJDeltaAlpha * getOsd_()->getEffectiveLambdaMatrix(getImpactBody())
             * getOsd_()->getDcJacobianInv(getImpactBody()))
      * tempReductionProjector * getOsd_()->getJacobian(getImpactBody());
  // ----------------------------------- Jacobian Tau
  Eigen::MatrixXd tempJDeltaTau;
  tempJDeltaTau.resize(getOsd_()->getDof(), 3);
  tempJDeltaTau.setZero();

  for(auto it = cache_.grfContainer.begin(); it != cache_.grfContainer.end(); ++it)
  {
    // auto it = cache_.grfContainer.find("l_sole");
    /*
        if(it->first == getImpactBody())
        {
          // Skip impact body
          continue;
        }
    */
    if(it->second.contact())
    {
      Eigen::MatrixXd temp;
      if(useLinearJacobian_())
      {
        temp.resize(3, 3);
      }
      else
      {
        temp.resize(6, 6);
      }
      temp.setZero();

      for(auto idx = cache_.grfContainer.begin(); idx != cache_.grfContainer.end(); ++idx)
      {
        temp += getOsd_()->getLambdaMatrix(it->first, idx->first)
                * getOsd_()->getLambdaMatrixInv(it->first, getImpactBody());
        // std::cout<<"temp is: "<<temp<<std::endl;
        // std::cout<<"temp is: "<<getOsd_()->getLambdaMatrix( it->first, idx->first) *
        // getOsd_()->getLambdaMatrixInv(idx->first, getImpactBody())<<std::endl;
      } // end of inner loop
      it->second.jacobianDeltaF = (1 / getImpactDuration()) * temp
                                  * getOsd_()->getEffectiveLambdaMatrix(getImpactBody())
                                  * getOsd_()->getDcJacobianInv(getImpactBody()) * tempReductionProjector
                                  * getOsd_()->getJacobian(getImpactBody());

      tempJDeltaTau += getOsd_()->getJacobian(it->first).transpose() * temp;
      // std::cout<<"tempJdeltaTau is: "<<tempJDeltaTau<<std::endl;
    }
  } // end of temp J delta tau

  /*
  cache_.jacobianDeltaTau = -(1 + getCoeRes_()) * (1 / getImpactDuration()) * tempJDeltaTau
                            * getOsd_()->getEffectiveLambdaMatrix(getImpactBody())
                            * getOsd_()->getDcJacobianInv(getImpactBody()) * tempProjector
                            * getOsd_()->getJacobian(getImpactBody());
  */

  cache_.jacobianDeltaTau =
      (getOsd_()->getJacobian(getImpactBody()).transpose() + (1 / getImpactDuration()) * tempJDeltaTau)
      * getOsd_()->getEffectiveLambdaMatrix(getImpactBody()) * getOsd_()->getDcJacobianInv(getImpactBody())
      * tempReductionProjector * getOsd_()->getJacobian(getImpactBody());
}
void mi_impactPredictor::printInfo() 
{

  std::cout<<"mi_impactPredictor with impactBody: "<<getImpactBody()<<", with OSD model: "<<std::endl;
  getOsd_()->printInfo();
	  

}
