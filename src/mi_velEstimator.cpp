/* Copyright 2019 CNRS-UM LIRMM
 *
 * \author Yuquan Wang, Arnaud Tanguy
 *
 *
 *
 * mc_impact_predictor is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * mc_impact_predictor is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mc_impact_predictor. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "mi_velEstimator.h"

namespace mc_impact
{

mi_velEstimator::mi_velEstimator(const mc_rbdyn::Robot & simRobot,
                               const std::shared_ptr<mi_osd> osdPtr,
			       const std::shared_ptr<mi_qpEstimator> qpEstimator,
                               const struct qpEstimatorParameter params)
: simRobot_(simRobot), osdPtr_(osdPtr), qpEstimator_(qpEstimator), params_(params), solverTime_(0), structTime_(0.0)
{

  // (0) Initialize the centroidal momentum calculator
  cmmPtr_ = std::make_shared<rbd::CentroidalMomentumMatrix>(osdPtr_->getRobot().mb());
  // (1) initilize the Impact models
  for(std::map<std::string, Eigen::Vector3d>::const_iterator idx = params.impactNameAndNormals.begin();
      idx != params.impactNameAndNormals.end(); ++idx)
  {
    ImpactModelParams params;
    params.coeF = params_.coeFrictionDeduction;
    params.iBodyName = idx->first;
    params.inertial_surfaceNormal = idx->second;
    params.coeR = params_.coeRes;
    params.timeStep = params_.timeStep;
    params.iDuration = params_.impactDuration;
    params.dim = params_.dim;
    params.useBodyJacobian = params_.impactModelBodyJacobian;

    impactModels_[idx->first] = std::make_shared<mc_impact::mi_impactModel>(getSimRobot(), params);
  }
  // Initialize the TwoDim frictional impact dynamics model
  TwoDimModelBridgeParams newTwoDimModelParams;
  newTwoDimModelParams.modelParams = getImpactModel("r_wrist")->getParams();

  newTwoDimModelParams.name = "PredictorTwoDimModelSimRobot";
  newTwoDimModelParams.useVirtualContact = false;
  newTwoDimModelParams.useComVel = false;
  twoDimFidModelPtr_ = std::make_shared<mc_impact::TwoDimModelBridge>(getSimRobot(), newTwoDimModelParams);

  // Try to use all the end-effectors from the OSD
  for(auto & ee : getOsd()->getEes())
  {
    bool isOsdEe = true;
    addEndeffector_(ee, isOsdEe);
  }

  if(params_.useJsd)
    eqConstraints_.emplace_back(
        std::make_shared<mc_impact::mi_velJsdEquality>(getOsd(), getQpEstimator_()));

  if(params_.useOsd)
    eqConstraints_.emplace_back(
        std::make_shared<mc_impact::mi_velInvOsdEquality>(getOsd(), getQpEstimator_(), static_cast<int>(getImpactModels().size())));

  if(params_.useImpulseBalance)
    eqConstraints_.emplace_back(std::make_shared<mc_impact::mi_velBalance>(getOsd(), getImpactModels(), cmmPtr_, getFidModel()));

  /*
  if(params_.useContactConstraint)
    eqConstraints_.emplace_back(std::make_shared<mc_impact::mi_contactConstraint>(getOsd(), getImpactModels(), endEffectors_));
    */

  for(std::map<std::string, Eigen::Vector3d>::const_iterator idx = params.impactNameAndNormals.begin();
      idx != params.impactNameAndNormals.end(); ++idx)
  {
    eqConstraints_.emplace_back(
        std::make_shared<mc_impact::mi_velIniEquality>(getOsd(), getImpactModel(idx->first), getEeNum()));
  }
  
  vector_A_dagger_.resize(params.impactNameAndNormals.size());

  std::cout << "Created QP estimator constraint. " << std::endl;

  initializeQP_();
  std::cout << "the QP-based impulse estimator is created. " << std::endl;

  }

mi_velEstimator::~mi_velEstimator()
{
  if(logEntries_.size() > 0)
  {
    assert(hostCtlPtr_ != nullptr);
    removeImpulseEstimations_();
  }

  if(guiEntries_.size() > 0)
  {
    assert(hostCtlPtr_ != nullptr);
    removeMcRtcGuiItems();
  }
}

void mi_velEstimator::initializeQP_()
{
  numVar_ = getDof();

  numEq_ = 0;
  for(auto idx = eqConstraints_.begin(); idx != eqConstraints_.end(); ++idx)
  {
    numEq_ += (*idx)->nrEq();
  }

  solver_.resize(getNumVar_(), getNumEq_(), Eigen::lssol::QP2);

  xl_.resize(getNumVar_());
  xu_.resize(getNumVar_());
  xl_ = xl_.setOnes() * -std::numeric_limits<double>::infinity();
  xu_ = xu_.setOnes() * std::numeric_limits<double>::infinity();

  p_.resize(getNumVar_());
  p_.setZero();

  Q_.resize(getNumVar_(), getNumVar_());
  Q_ = Q_.setIdentity() * getQweight();

  C_.resize(getNumEq_(), getNumVar_());
  C_.setZero();
  cu_.resize(getNumEq_());
  cu_.setZero();
  cl_.resize(getNumEq_());
  cl_.setZero();

  jointVelJump_.resize(getDof());
  tauJump_.resize(getDof());
  jointVelJump_.setZero();
  tauJump_.setZero();

  jacobianDeltaAlpha_.resize(getDof(), getDof());
  jacobianDeltaAlpha_.setZero();

  jacobianTwoDeltaAlpha_.resize(getDof(), getDof());
  jacobianTwoDeltaAlpha_.setZero();

  jacobianDeltaTau_.resize(getDof(), getDof());
  jacobianDeltaTau_.setZero();

  // Resize each element of vector_A_dagger_, where each of them correspond to one impact
  /*
  for(auto & A_dagger : vector_A_dagger_)
  {
    A_dagger.resize(getDof() + getEstimatorParams().dim * getEeNum(), 3);

    A_dagger.setZero();
  }
  */
  std::cout << "Reset LSSOL QP estimator variables. " << std::endl;
}

void mi_velEstimator::solveWeightedEqQp_(const Eigen::MatrixXd & Q_,
                                        const Eigen::VectorXd & p_,
                                        const Eigen::MatrixXd & C_,
                                        const Eigen::VectorXd & cu_,
                                        Eigen::VectorXd & solution)
{
  Eigen::MatrixXd kkt;
  int kktDim = getNumVar_() + static_cast<int>(C_.rows());

  kkt.resize(kktDim, kktDim);
  kkt.setZero();

  kkt.block(0, 0, getNumVar_(), getNumVar_()) = Q_;
  // kkt.block(0, 0, getNumVar_(), getNumVar_()).setIdentity(getNumVar_(), getNumVar_());

  // kkt.block(0, 0, getNumVar_(), getNumVar_())=  kkt.block(0, 0, getNumVar_(), getNumVar_())*2;

  kkt.block(getNumVar_(), 0, C_.rows(), C_.cols()) = C_;
  kkt.block(0, getNumVar_(), C_.cols(), C_.rows()) = C_.transpose();

  Eigen::VectorXd b;
  b.resize(kktDim);
  b.setZero();
  b.segment(getNumVar_(), cu_.rows()) = cu_;

  // solve the least squares problem:
  // Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_kkt(kkt);
  // solution.resize(kktDim);
  // solution = lu_decomp_kkt.inverse()*cu_;

  // Eigen::MatrixXd tempInverse =  kkt.completeOrthogonalDecomposition().pseudoInverse();

  // Suppose kkt has full column rank, we perform cholesky decomposition to compute pinv = (A^*A)^{-1}A^*
  Eigen::MatrixXd tempInverse = (kkt.transpose() * kkt).inverse() * kkt.transpose();

  // We discard the Lagrange multipliers associated with C_;
  solution = (tempInverse * b).segment(0, getNumVar_());
  // solution = kkt.completeOrthogonalDecomposition().pseudoInverse()*b;

  // std::cout<<"tempInverse size is: "<<tempInverse.rows() <<", "<<tempInverse.cols()<<std::endl;

  // std::cout<<"test 2" <<std::endl;
  for(size_t idx = 0; idx < impactModels_.size(); ++idx)
  {
    vector_A_dagger_[idx] =
        tempInverse.block(0, tempInverse.cols() - 3 * static_cast<int>(idx + 1), vector_A_dagger_[idx].rows(), 3);
  }

  // std::cout<<"test 3" <<std::endl;
  // A_dagger_ =  //tempInv_= tempInverse.block(0, tempInverse.cols() - 3, tempInverse.rows(), 3);
  tempInv_ = tempInverse;
}

void mi_velEstimator::solveEqQp_(const Eigen::MatrixXd & Q_,
                                const Eigen::VectorXd & p_,
                                const Eigen::MatrixXd & C_,
                                const Eigen::VectorXd & cu_,
                                Eigen::VectorXd & solution)
{
  Eigen::MatrixXd kkt;
  int kktDim = getNumVar_() + static_cast<int>(C_.rows());

  kkt.resize(kktDim, kktDim);
  kkt.setZero();

  // kkt.block(0, 0, getNumVar_(), getNumVar_()) = Q_;
  kkt.block(0, 0, getNumVar_(), getNumVar_()).setIdentity(getNumVar_(), getNumVar_());

  // kkt.block(0, 0, getNumVar_(), getNumVar_())=  kkt.block(0, 0, getNumVar_(), getNumVar_())*2;

  kkt.block(getNumVar_(), 0, C_.rows(), C_.cols()) = C_;
  kkt.block(0, getNumVar_(), C_.cols(), C_.rows()) = C_.transpose();

  Eigen::VectorXd b;
  b.resize(kktDim);
  b.setZero();
  b.segment(getNumVar_(), cu_.rows()) = cu_;

  // solve the least squares problem:
  // Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_kkt(kkt);
  // solution.resize(kktDim);
  // solution = lu_decomp_kkt.inverse()*cu_;

  // Eigen::MatrixXd tempInverse =  kkt.completeOrthogonalDecomposition().pseudoInverse();

  // Suppose kkt has full column rank, we perform cholesky decomposition to compute pinv = (A^*A)^{-1}A^*
  Eigen::MatrixXd tempInverse = (kkt.transpose() * kkt).inverse() * kkt.transpose();

  solution = tempInverse * b;
  // solution = kkt.completeOrthogonalDecomposition().pseudoInverse()*b;

  // std::cout<<"tempInverse size is: "<<tempInverse.rows() <<", "<<tempInverse.cols()<<std::endl;

  // std::cout<<"test 2" <<std::endl;
  for(size_t idx = 0; idx < impactModels_.size(); ++idx)
  {
    vector_A_dagger_[idx] =
        tempInverse.block(0, tempInverse.cols() - 3 * static_cast<int>(idx + 1), vector_A_dagger_[idx].rows(), 3);
  }

  // std::cout<<"test 3" <<std::endl;
  // A_dagger_ =  //tempInv_= tempInverse.block(0, tempInverse.cols() - 3, tempInverse.rows(), 3);
  tempInv_ = tempInverse;
}
void mi_velEstimator::updateImpactModels_()
{
  for(auto idx = getImpactModels().begin(); idx != getImpactModels().end(); ++idx)
  {
    //(getImpactModel(const_cast<std::string & >(idx->first)))->update(idx->second);
    //(getImpactModel(idx->first))->update();
    idx->second->update();
  }

  if(getOsd()->useBodyJacobian())
  {
    // The two dim model by defaults computes in the inertial frame
    // Convert the body frame variables to the inertial frame
    const Eigen::Matrix3d & rotation = getSimRobot().bodyPosW("r_wrist").rotation();
    twoDimFidModelPtr_->update(rotation * getImpactModel("r_wrist")->getSurfaceNormal(),
                               rotation * getImpactModel("r_wrist")->getEeVelocity());
  }
  else
  {
    // Inertial Frame
    twoDimFidModelPtr_->update(getImpactModel("r_wrist")->getSurfaceNormal(),
                               getImpactModel("r_wrist")->getEeVelocity());
  }
}

void mi_velEstimator::updateImpactModels_(const std::map<std::string, Eigen::Vector3d> & surfaceNormals)
{
  for(std::map<std::string, Eigen::Vector3d>::const_iterator idx = surfaceNormals.begin(); idx != surfaceNormals.end();
      ++idx)
  {
    //(getImpactModel(const_cast<std::string & >(idx->first)))->update(idx->second);
    (getImpactModel(idx->first))->update(idx->second);
  }
}
void mi_velEstimator::update()
{
  updateImpactModels_();
  update_();
}

void mi_velEstimator::update(const std::map<std::string, Eigen::Vector3d> & surfaceNormals)
{
  if(surfaceNormals.size() != impactModels_.size())
  {
    throw std::runtime_error(std::string("mi_velEstimator-update: surfaceNormals size(")
                             + std::to_string(static_cast<int>(surfaceNormals.size()))
                             + std::string(") does not match impact predictor impact number (")
                             + std::to_string(impactModels_.size()) + std::string(")."));
  }

  updateImpactModels_(surfaceNormals);
  update_();
}

void mi_velEstimator::update_()
{

  // Update the CMM matrix
  cmmPtr_->sComputeMatrix(osdPtr_->getRobot().mb(), osdPtr_->getRobot().mbc(), osdPtr_->getRobot().com());

  auto startStruct = std::chrono::high_resolution_clock::now();
  int count = 0;
  // Update the constraints
  for(auto & eq : eqConstraints_)
  {
    eq->update();

    // Range space
    C_.block(count, 0, eq->nrEq(), getNumVar_()) = eq->AEq();

    // Lower bounds: lb == ub for equalities
    cl_.segment(count, eq->nrEq()) = eq->bEq();

    // Upper bounds:
    cu_.segment(count, eq->nrEq()) = eq->bEq();

    count += eq->nrEq();
  }

  Eigen::VectorXd solutionVariables;
  solutionVariables.resize(getNumVar_());
  solutionVariables.setZero();

  auto stopStruct = std::chrono::high_resolution_clock::now();

  auto durationStruct = std::chrono::duration_cast<std::chrono::microseconds>(stopStruct - startStruct);

  structTime_ = static_cast<double>(durationStruct.count());

  // Update Q with the current mass matrix
  auto startSolve = std::chrono::high_resolution_clock::now();
  if(getEstimatorParams().useLagrangeMultiplier)
  {

    if(getEstimatorParams().testWeightedQp)
    {
      solveWeightedEqQp_(Q_, p_, C_, cu_, solutionVariables);
    }
    else
    {
      solveEqQp_(Q_, p_, C_, cu_, solutionVariables);
    }
  }
  else
  {
    solver_.solve(xl_, xu_, Q_, p_, C_, cl_, cu_);
    solutionVariables = solver_.result();
  }

  auto stopSolve = std::chrono::high_resolution_clock::now();
  auto durationSolve = std::chrono::duration_cast<std::chrono::microseconds>(stopSolve - startSolve);

  solverTime_ = static_cast<double>(durationSolve.count());

  // std::cout<<"test 1" <<std::endl;
  // Eigen::MatrixXd tempJac = getImpactModel()->getProjector();
  jointVelJump_ = solutionVariables;

  double mass_inv = 1.0 / getSimRobot().mass();

  const Eigen::MatrixXd & cmmMatrix = getCmm()->matrix();

  // First three entries correspond to the angular part
  comVelJump_ = mass_inv * cmmMatrix.block(3, 0, 3, getDof()) * jointVelJump_;

  tauJump_.setZero();

  jacobianDeltaTau_.resize(getDof(), getDof());

  jacobianDeltaTau_.setZero();

  readEeJacobiansSolution_(solutionVariables);

  if(getEstimatorParams().useLagrangeMultiplier)
  {
    unsigned int iiA = 0;
    for(auto impactIdx = impactModels_.begin(); impactIdx != impactModels_.end(); ++impactIdx, ++iiA)
    {
      jacobianDeltaAlpha_ = vector_A_dagger_[iiA].block(0, 0, getDof(), 3) * impactIdx->second->getProjector();
      jacobianTwoDeltaAlpha_ = vector_A_dagger_[iiA].block(0, 0, getDof(), 3) * impactIdx->second->getProjectorTwo();
    }
  }

}

const endEffector & mi_velEstimator::getEndeffector(const std::string & name)
{
  return getEndeffector_(name);
}

bool mi_velEstimator::osdContactEe_(const std::string & eeName)
{
  return getOsd()->hasEndeffector(eeName);
}

endEffector & mi_velEstimator::getEndeffector_(const std::string & name)
{
  auto opt = endEffectors_.find(name);
  if(opt != (endEffectors_.end()))
  {
    return opt->second;
  }
  else
  {
    throw std::runtime_error(std::string("getEndeffector: '") + name + std::string("' is not found."));
  }
}

const std::shared_ptr<mc_impact::mi_impactModel> mi_velEstimator::getImpactModel(const std::string & eeName)
{
  auto idx = impactModels_.find(eeName);
  if(idx != (impactModels_.end()))
  {
    return idx->second;
  }
  else
  {
    throw std::runtime_error(std::string("getImpactModel: '") + eeName + std::string("' is not found."));
  }
}

int mi_velEstimator::nameToIndex_(const std::string & eeName)
{
  auto tempEe = endEffectors_.find(eeName);
  if(tempEe != endEffectors_.end())
    return tempEe->second.uniqueIndex;
  else
  {
    std::string error_msg = std::string("qpEstimator::nameToIndex_: ee-") + eeName + std::string(": does not exist.");
    throw std::runtime_error(error_msg);
  }
}

bool mi_velEstimator::addEndeffector_(const std::string & eeName, const bool & fromOsdModel)
{

  auto opt = endEffectors_.find(eeName);
  if(opt != (endEffectors_.end()))
  {
    std::cout << "Endeffector " << eeName << " already exists." << std::endl;
    return false;
  }

  // addOptVariables_(eeName, getDim());
  Eigen::Vector3d tempForce = Eigen::Vector3d::Zero();

  Eigen::MatrixXd tempJ;
  tempJ.resize(getEstimatorParams().dim, getDof());
  tempJ.setZero();

  sva::ForceVecd tempWrench = sva::ForceVecd::Zero();

  int eeIndex = 0;

  if(fromOsdModel)
    eeIndex = getOsd()->nameToIndex_(eeName);
  else
    eeIndex = static_cast<int>(endEffectors_.size());

  endEffectors_[eeName] = {eeIndex, tempForce, tempForce, tempForce, tempForce, tempForce, tempJ, tempWrench};

  std::cout << "Qp Estimator: Adding end-effector: " << eeName << ", with index: " << eeIndex << std::endl;
  return true;
}

void mi_velEstimator::print(const std::string & eeName)
{
  auto & tempEe = getEndeffector(eeName);
  std::cout << "Endeffector " << eeName << " predicted vel jump: " << tempEe.eeVJump.transpose()
            << ", predicted impulse force: " << tempEe.estimatedAverageImpulsiveForce.transpose()
            << ", predicted impulse: " << tempEe.estimatedImpulse.transpose() << std::endl;
}
void mi_velEstimator::print() const
{
  /*
  std::cout<<"The QP estimator params are: "<<std::endl<<"Dim: " <<getImpactModel()->getDim() <<", Dof: "<<getDof()<<",
  coeR: "<<getImpactModel()->getCoeRes()<<", coeF: "<<getImpactModel()->getCoeFricDe()<<", impact duration:
  "<<getImpactModel()->getImpactDuration()<<". "<<std::endl;
*/
  std::cout << red << "The QP estimator: " << getEstimatorParams().name
            << " has an OSD model with the end-effectors: " << cyan;
  for(auto idx = getOsd()->getEes().begin(); idx != getOsd()->getEes().end(); ++idx)
  {
    std::cout << *idx << " ";
  }
  std::cout << reset << std::endl;

  std::cout << red << "The QP estimator: " << getEstimatorParams().name
            << " has an  OSD model  with established contacts: " << green;
  for(auto idx = getOsd()->getContactEes().begin(); idx != getOsd()->getContactEes().end(); ++idx)
  {
    std::cout << *idx << " ";
  }
  std::cout << reset << std::endl;
}



void mi_velEstimator::readEeJacobiansSolution_(const Eigen::VectorXd & solutionVariables)
{

  // We fill the following for each end-effector:
  // (a) force Jump Jacobian
  // (b) impulse
  // (c) force Jump
  //
  // We also fill:
  // Joint torque Jump

  for(auto idx = endEffectors_.begin(); idx != endEffectors_.end(); ++idx)
  {
    int eeIndex = nameToIndex_(idx->first);
    int location = getEstimatorParams().dim * eeIndex;

    // auto & tempEe =  getEndeffector_(*idx);
    double inv_dt = 1.0 / (getImpactModels().begin()->second->getParams().iDuration);

    //idx->second.estimatedImpulse = solutionVariables.segment(getDof() + location, getEstimatorParams().dim);

    //idx->second.estimatedAverageImpulsiveForce = idx->second.estimatedImpulse * inv_dt;
    idx->second.rssForce = inv_dt * getOsd()->getEquivalentMass(idx->first) * getOsd()->getJacobian(idx->first) * getJointVelJump();  
    Eigen::MatrixXd tempJ;

    // The Jacobian may exist in two places
    if(osdContactEe_(idx->first))
      tempJ = getOsd()->getJacobian(idx->first);
    else
      tempJ = getImpactModel(idx->first)->getJacobian();

    idx->second.eeVJump = tempJ * jointVelJump_;

    // std::cout<<"A_dagger_ size is: "<<A_dagger_.rows() <<", "<<A_dagger_.cols()<<std::endl;
    // std::cout<<"Dof is: "<<getDof()<<std::endl;
    // std::cout<<"location is: "<<location<<std::endl;

    // std::cout<<"tempA_dagger_ee is: "<<std::endl<<tempA_dagger_ee<<std::endl;
    if(getEstimatorParams().useLagrangeMultiplier)
    {

      // int kktDim =  getNumVar_() + static_cast<int>(C_.rows());
      // Eigen::VectorXd b;
      // b.resize(kktDim); b.setZero();
      // b.segment(getNumVar_(), cu_.rows() ) = cu_;

      // auto solution = tempInv_*b;

      // double inv_dt = 1.0/getImpactModel()->getImpactDuration();
      // Eigen::MatrixXd test_A = tempInv_.block(getDof() + location, 0, 3, 3)  ;
      // tempEe.checkForce = test_A*getImpactModel()->getEeVelocityJump()*inv_dt;
      // tempEe.checkForce = inv_dt*solution.segment(getDof() + location, 3);

      idx->second.jacobianDeltaF.resize(3, getDof());
      idx->second.jacobianDeltaF.setZero();

      // for (int i?i = 0; ii< static_cast<int>(impactModels_.size()); ++ii )
      unsigned int iiA = 0;
      Eigen::MatrixXd tempJ_T = tempJ.transpose();
      for(auto impactIdx = impactModels_.begin(); impactIdx != impactModels_.end(); ++impactIdx, ++iiA)
      {

        Eigen::Matrix3d tempA_dagger_ee = vector_A_dagger_[iiA].block(getDof() + location, 0, 3, 3);

        idx->second.jacobianDeltaF += tempA_dagger_ee * impactIdx->second->getProjector();

        idx->second.checkForce += inv_dt * idx->second.jacobianDeltaF * impactIdx->second->getJointVel();

        jacobianDeltaTau_ += tempJ_T * tempA_dagger_ee * impactIdx->second->getProjector();
      }

      tauJump_ += tempJ_T * idx->second.estimatedAverageImpulsiveForce;

      // tempEe.checkForce = inv_dt* tempInv_.block(getDof() + location, tempInv_.cols() - 3, 3,
      // 3)*getImpactModel()->getEeVelocityJump();

      // tempEe.checkForce = inv_dt* tempInv_.block(getDof() + location, tempInv_.cols() , 3, tempInv_.cols())*b;
      // tempEe.checkForce(0) = (b.segment(0, b.rows()-3)).norm();
      // std::cout<<"b is: " <<std::endl<<b.transpose()<<std::endl;
      /*
            tempEe.checkForce(0) = (
                //jointVelJump_  -  tempInv_.block(0, tempInv_.cols() - 3, getDof(), 3)*
         getImpactModel()->getEeVelocityJump()
                //jointVelJump_  -  tempInv_.block(0, tempInv_.cols() - 3, getDof(), 3)* tempJac*
         getImpactModel()->getJointVel() jointVelJump_  -  A_dagger_.block(0, 0, getDof(), 3)* tempJac*
         getImpactModel()->getJointVel()
                ).norm();
      */

      // std::cout<<"the impulse difference is: "<<(tempEe.estimatedImpulse -
      // tempA_dagger_ee*getImpactModel()->getEeVelocityJump()).norm()<<std::endl;
      // Transform of the endEffector.

      auto X_0_c = getSimRobot().bodyPosW(idx->first);
      auto translation = X_0_c.translation() - getSimRobot().com();
      // R_0_pi
      auto rotation = X_0_c.rotation();

      Eigen::Matrix3d torqueMatrix = getOsd()->crossMatrix(translation);
      Eigen::Matrix3d forceMatrix = Eigen::Matrix3d::Identity();

      if(getOsd()->useBodyJacobian())
      {
        torqueMatrix = torqueMatrix * rotation;
        forceMatrix = rotation;
      }


    } // end of Lagrange Multipliers

    // tempEe.jacobianDeltaF = tempA_dagger_ee;

    // std::cout<<"test "<<std::endl;
    // std::cout<<"test "<<std::endl;

    // std::cout<<"tempJ_T size is: "<<tempJ_T.rows() <<", "<<tempJ_T.cols()<<std::endl;

    // std::cout<<"test "<<std::endl;

    // std::cout<<"test "<<std::endl;

  } // end of the for end-effector loop.
}

void mi_velEstimator::logImpulseEstimations()
{

  const std::string & qpName = getEstimatorParams().name;
  logEntries_.emplace_back(qpName + "_" + "JointVelJump");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getJointVelJump(); });

  logEntries_.emplace_back(qpName + "_" + "COMVelJump");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getCOMVelJump(); });

  logEntries_.emplace_back(qpName + "_" + "CentroidalMomentum_simRobot");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                      [this]() { return getOsd()->getSimulatedCentroidalMomentum(); });

  logEntries_.emplace_back(qpName + "_" + "CentroidalMomentumD_simRobot");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                      [this]() { return getOsd()->getSimulatedCentroidalMomentumD(); });

  logEntries_.emplace_back(qpName + "_" + "CentroidalMomentumJump_FID_linear");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]()->Eigen::Vector3d { return getOsd()->getRobot().mass() * getFidModel()->getRobotPostImpactStates().linearVelJump; });

  logEntries_.emplace_back(qpName + "_" + "CentroidalMomentumJump_FID_angular");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]()->Eigen::Vector3d { return getFidModel()->getRobotCentroidalInertia() * getFidModel()->getRobotPostImpactStates().anguleVelJump;
 });

  logEntries_.emplace_back(qpName + "_" + "JointTorqueJump");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getTauJump(); });

  // (1) Loop over the end-effectors:

  for(auto & ee : getOsd()->getEes())
  {
    // (1.1) Force jump
    logEntries_.emplace_back(qpName + "_" + ee + "_" + "ForceJump");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this, ee]() -> Eigen::Vector3d {
      return getEndeffector(ee).estimatedAverageImpulsiveForce;
    });

    // Impulse
    logEntries_.emplace_back(qpName + "_" + ee + "_" + "Impulse");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this, ee]() -> Eigen::Vector3d {
      return getEndeffector(ee).estimatedImpulse;
    });


    // This one should be the same with 'ForceJump'
    logEntries_.emplace_back(qpName + "_" + ee + "_" + "ForceJump_debug");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                        [this, ee]() -> Eigen::Vector3d { return getEndeffector(ee).checkForce; });

    logEntries_.emplace_back(qpName + "_" + ee + "_" + "ForceJump_debug_rss");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                        [this, ee]() -> Eigen::Vector3d { return getEndeffector(ee).rssForce; });


    // (1.2) Ee velocity jump
    logEntries_.emplace_back(qpName + "_" + ee + "_" + "eeVelJump");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                        [this, ee]() -> Eigen::Vector3d { return getEndeffector(ee).eeVJump; });
  }

  // (2) Loop over the impact bodies:
  for(auto & eePair : getImpactModels())
  {
    const std::string & eeName = eePair.first;
    // (2.1) Estimated end-effector induced impulse joint torque
    logEntries_.emplace_back(qpName + "_" + eeName + "_" + "jointTorqueJump");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this, eeName]() -> Eigen::Vector3d {
      return getImpactModel(eeName)->getJacobian().transpose() * getHostCtl_()->robot().bodyWrench(eeName).force();
    });

    // (2.2) Predicted ee velocity jump (using coefficient of restitution)
    logEntries_.emplace_back(qpName + "_" + eeName + "_" + "ImpactModel_eeVelJump");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                        [this, eeName]() { return getImpactModel(eeName)->getEeVelocityJump(); });

    // (2.3) ee velocity
    logEntries_.emplace_back(qpName + "_" + eeName + "_" + "ImpactModel_eeVel");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                        [this, eeName]() { return getImpactModel(eeName)->getEeVelocity(); });

    // (2.4) impact normal orientation
    logEntries_.emplace_back(qpName + "_" + eeName + "_" + "ImpactModel_impactNormal");
    getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                        [this, eeName]() { return getImpactModel(eeName)->getSurfaceNormal(); });
  }

  // Log the frictional impact dynamics entries.
  getFidModel()->logImpulseEstimations();
}

void mi_velEstimator::removeImpulseEstimations_()
{
  assert(getHostCtl_() != nullptr);

  getFidModel()->removeImpulseEstimations();

  for(auto & name : logEntries_)
  {
    getHostCtl_()->logger().removeLogEntry(name);
  }
}

void mi_velEstimator::addMcRtcGuiItems()
{

  mc_rtc::gui::ArrowConfig surfaceXConfig({0., (153.0 / 255.0), (153.0 / 255.0)});
  surfaceXConfig.start_point_scale = 0.0;
  surfaceXConfig.end_point_scale = 0.0;

  double arrowLengthScale = 0.01;

  // Loop over all the end-effectors:

  for(auto & eeName : getOsd()->getEes())
  {

    guiEntries_.emplace_back(eeName + "_ForceJump");

    getHostCtl_()->gui()->addElement(
        {getEstimatorParams().name},
        mc_rtc::gui::Arrow(guiEntries_.back(), surfaceXConfig,
                           [this, eeName]() -> Eigen::Vector3d {
                             // start of the arrow
                             auto X_0_s = getHostCtl_()->robot().bodyPosW(eeName);
                             return X_0_s.translation();
                           },
                           [this, eeName, arrowLengthScale]() -> Eigen::Vector3d {
                             // End of the arrow
                             auto X_0_s = getHostCtl_()->robot().bodyPosW(eeName);
                             // Note that (1) the forceJump is defined in the local frame.
                             // (2) the rotaiton =  X_0_s.rotation().transpose();
                             // (3) We are plotting the forceJump in the inertial frame.

                             Eigen::Vector3d eeForceJump;

                             if(getOsd()->useBodyJacobian())
                             {
                               // Force are computed in the body frame.
                               eeForceJump = X_0_s.translation()
                                             + arrowLengthScale * X_0_s.rotation().transpose()
                                                   * getEndeffector(eeName).estimatedAverageImpulsiveForce;
                             }
                             else
                             {
                               // Force are computed in the inertial frame
                               eeForceJump = X_0_s.translation()
                                             + arrowLengthScale * getEndeffector(eeName).estimatedAverageImpulsiveForce;
                             }
                             return eeForceJump;
                           }));
  }
}

void mi_velEstimator::removeMcRtcGuiItems()
{
  assert(getHostCtl_() != nullptr);

  for(auto & name : guiEntries_)
  {
    getHostCtl_()->gui()->removeElement({getEstimatorParams().name}, name);
  }
}

} // namespace mc_impact