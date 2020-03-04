#include "mi_qpEstimator.h"

namespace mc_impact
{

mi_qpEstimator::mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
                               const std::shared_ptr<mi_osd> osdPtr,
                               const struct qpEstimatorParameter params)
: simRobot_(simRobot), osdPtr_(osdPtr), params_(params), solverTime_(0), structTime_(0.0)
{
	
  // (1) initilize the Impact models
  for(std::map<std::string, Eigen::Vector3d>::const_iterator idx = params.impactNameAndNormals.begin();
      idx != params.impactNameAndNormals.end(); ++idx)
  {
    impactModels_[idx->first] = std::make_shared<mc_impact::mi_impactModel>(
        getSimRobot(), idx->first, idx->second, params_.impactDuration, params_.timeStep, params_.coeFrictionDeduction,
        params_.coeRes, params_.dim);
  }
  // (2) Add the end-effectors: first OSD endeffectors, then the impact model endeffectors
  /*
  for(auto idx = getOsd()->getContactEes().begin(); idx != getOsd()->getContactEes().end(); ++idx)
  {
    bool isOsdEe = true;
    addEndeffector_(*idx, isOsdEe);
  }
  for(auto idx = getImpactModels().begin(); idx != getImpactModels().end(); ++idx)
  {

    bool isOsdEe = false;
    addEndeffector_(idx->first, isOsdEe);
  }

  */
 
  // Try to use all the end-effectors from the OSD
  for(auto & ee: getOsd()->getEes())
  {
    bool isOsdEe = true;
    addEndeffector_(ee, isOsdEe);
  }

  if(params_.useJsd)
    eqConstraints_.emplace_back(std::make_shared<mc_impact::mi_jsdEquality>(getOsd(), getImpactModels(), endEffectors_));

  if(params_.useOsd)
    eqConstraints_.emplace_back(
        std::make_shared<mc_impact::mi_invOsdEquality>(getOsd(), static_cast<int>(getImpactModels().size())));

  if(params_.useImpulseBalance)
    eqConstraints_.emplace_back(
        std::make_shared<mc_impact::mi_balance>(getOsd(), getImpactModels()));


  for(std::map<std::string, Eigen::Vector3d>::const_iterator idx = params.impactNameAndNormals.begin();
      idx != params.impactNameAndNormals.end(); ++idx)
  {
    // eqConstraints_.push_back(std::make_shared<mi_iniEquality>(getOsd(),
    // getImpactModel(const_cast<std::string&>(*idx)).get(), false));
    eqConstraints_.emplace_back(
        std::make_shared<mc_impact::mi_iniEquality>(getOsd(), getImpactModel(idx->first), getEeNum()));
  }
  /*
    for (std::vector<std::string>::const_iterator idx = params.impactBodyNames.begin();
    idx!=params.impactBodyNames.end(); ++idx)
    {
     impactModels_[*idx] = std::make_shared<mi_impactModel>(getSimRobot(), getOsd(), *idx, params_.impactDuration,
    params_.timeStep, params_.coeFrictionDeduction, params_.coeRes, params_.dim);

    //eqConstraints_.push_back(std::make_shared<mi_iniEquality>(getOsd(),
    getImpactModel(const_cast<std::string&>(*idx)).get(), false));
    eqConstraints_.push_back(std::make_shared<mi_iniEquality>(getOsd(), getImpactModel(*idx).get(), false));
    }
  */
  vector_A_dagger_.resize(params.impactNameAndNormals.size());

  std::cout << "Created QP estimator constraint. " << std::endl;

  initializeQP_();
  std::cout << "the QP-based impulse estimator is created. " << std::endl;
}

void mi_qpEstimator::initializeQP_()
{
  numVar_ = getDof() + 3 * getEeNum();

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
  jacobianDeltaTau_.resize(getDof(), getDof());
  jacobianDeltaTau_.setZero();

  // Resize each element of vector_A_dagger_, where each of them correspond to one impact
  for(auto & A_dagger : vector_A_dagger_ )
  {
    A_dagger.resize(
        getDof() + getEstimatorParams().dim * getEeNum(), 3);

    A_dagger.setZero();
  }
  std::cout << "Reset LSSOL QP estimator variables. " << std::endl;
}

void mi_qpEstimator::solveEqQp_(const Eigen::MatrixXd & Q_,
                                const Eigen::VectorXd & p_,
                                const Eigen::MatrixXd & C_,
                                const Eigen::VectorXd & cu_,
                                Eigen::VectorXd & solution)
{
  Eigen::MatrixXd kkt;
  int kktDim = getNumVar_() + static_cast<int>(C_.rows());

  kkt.resize(kktDim, kktDim);
  kkt.setZero();

  //kkt.block(0, 0, getNumVar_(), getNumVar_()) = Q_;
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
void mi_qpEstimator::updateImpactModels_()
{
  for(auto idx = getImpactModels().begin(); idx != getImpactModels().end(); ++idx)
  {
    //(getImpactModel(const_cast<std::string & >(idx->first)))->update(idx->second);
    //(getImpactModel(idx->first))->update();
    idx->second->update();
  }
}



void mi_qpEstimator::updateImpactModels_(const std::map<std::string, Eigen::Vector3d> & surfaceNormals)
{
  for(std::map<std::string, Eigen::Vector3d>::const_iterator idx = surfaceNormals.begin(); idx != surfaceNormals.end();
      ++idx)
  {
    //(getImpactModel(const_cast<std::string & >(idx->first)))->update(idx->second);
    (getImpactModel(idx->first))->update(idx->second);
  }
}
void mi_qpEstimator::update()
{
  updateImpactModels_();
  update_();
}

void mi_qpEstimator::update(const std::map<std::string, Eigen::Vector3d> & surfaceNormals)
{
  if(surfaceNormals.size() != impactModels_.size())
  {
    throw std::runtime_error(std::string("mi_qpEstimator-update: surfaceNormals size(")
                             + std::to_string(static_cast<int>(surfaceNormals.size()))
                             + std::string(") does not match impact predictor impact number (")
                             + std::to_string(impactModels_.size()) + std::string(")."));
  }

  updateImpactModels_(surfaceNormals);
  update_();
}
void mi_qpEstimator::update_()
{

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
  
  Q_.block(0 , 0, getDof(), getDof()) = getOsd()->getMassMatrix();

  auto startSolve = std::chrono::high_resolution_clock::now();
  if(params_.useLagrangeMultiplier)
  {

    // Eigen::VectorXd b_temp = cu_.segment(getNumVar_(), C_.rows());
    // Eigen::MatrixXd A_temp = C_;
    solveEqQp_(Q_, p_, C_, cu_, solutionVariables);
    // solutionVariables = A_temp.completeOrthogonalDecomposition().pseudoInverse()*b_temp;
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
  jointVelJump_ = solutionVariables.segment(0, getDof());
  tauJump_.setZero();

  jacobianDeltaTau_.resize(getDof(), getDof());

  // std::cout<<"jacobianDeltaTau_ size is: "<<jacobianDeltaTau_.rows() <<", "<<jacobianDeltaTau_.cols()<<std::endl;
  jacobianDeltaTau_.setZero();

  readEeJacobiansSolution_(solutionVariables);
  
  if(getEstimatorParams().useLagrangeMultiplier)
  {
    unsigned int iiA = 0;
    for(auto impactIdx = impactModels_.begin(); impactIdx != impactModels_.end(); ++impactIdx, ++iiA)
    {
      jacobianDeltaAlpha_ = vector_A_dagger_[iiA].block(0, 0, getDof(), 3) * impactIdx->second->getProjector();
    }
  }

  /*
    if(getEstimatorParams().useLagrangeMultiplier)
      //jacobianDeltaAlpha_ = A_dagger_.block(0, 0, getDof(), 3)*tempJac;
      jacobianDeltaAlpha_ = tempInv_.block(0, tempInv_.cols() - 3, getDof(), 3)*tempJac;
  */
  // jacobianDeltaAlpha_ = A_dagger_.block(0, 0, getDof(), 3);
  // jacobianDeltaTau_ = jacobianDeltaTau_*tempJac;
/* 
  for(auto idx = endEffectors_.begin(); idx != endEffectors_.end(); ++idx)
  {
    idx->second.perturbedWrench.vector().setZero();

    for(auto ii = endEffectors_.begin(); ii != endEffectors_.end(); ++ii)
    {
      if(ii==idx)
  continue;
      sva::PTransformd X_0_ii = getSimRobot().bodyPosW(ii->first);
      sva::PTransformd X_0_idx = getSimRobot().bodyPosW(idx->first);
      sva::PTransformd X_ii_idx = X_0_idx*X_0_ii.inv();
      sva::ForceVecd tempWrench;
      tempWrench.force().setZero();
      tempWrench.force() = ii->second.estimatedAverageImpulsiveForce;
      tempWrench.couple().setZero();

      idx->second.perturbedWrench +=  X_ii_idx.dualMul(tempWrench);

    }


  }

 */ 

  calcPerturbedWrench_();

}

const endEffector & mi_qpEstimator::getEndeffector(const std::string & name)
{
  return getEndeffector_(name);
}

bool mi_qpEstimator::osdContactEe_(const std::string & eeName)
{
  return getOsd()->hasEndeffector(eeName);
}

endEffector & mi_qpEstimator::getEndeffector_(const std::string & name)
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

const std::shared_ptr<mc_impact::mi_impactModel> & mi_qpEstimator::getImpactModel(const std::string & eeName)
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

int mi_qpEstimator::nameToIndex_(const std::string & eeName)
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

bool mi_qpEstimator::addEndeffector_(const std::string & eeName, const bool & fromOsdModel)
{

  auto opt = endEffectors_.find(eeName);
  if(opt != (endEffectors_.end()))
  {
    std::cout << "Endeffector " << eeName << " already exists." << std::endl;
    return false;
  }

  // addOptVariables_(eeName, getDim());
  Eigen::Vector3d tempForce;
  // tempForce.resize(getImpactModel()->getDim());
  tempForce.setZero();

  Eigen::MatrixXd tempJ;
  tempJ.resize(getEstimatorParams().dim, getDof());
  tempJ.setZero();

  
  sva::ForceVecd tempWrench;
  tempWrench.force().setZero();
  tempWrench.couple().setZero();
  
  int eeIndex = 0;

  if(fromOsdModel)
    eeIndex = getOsd()->nameToIndex_(eeName);
  else
    eeIndex = static_cast<int>(endEffectors_.size());

  // optVariables_[name] = {dim, optVariables_.size() };
  endEffectors_[eeName] = {eeIndex, tempForce, tempForce, tempForce, tempForce, tempJ, tempWrench};

  std::cout << "Qp Estimator: Adding end-effector: " << eeName << ", with index: " << eeIndex << std::endl;
  return true;
}

void mi_qpEstimator::print(const std::string & eeName)
{
  auto & tempEe = getEndeffector(eeName);
  std::cout << "Endeffector " << eeName << " predicted vel jump: " << tempEe.eeVJump.transpose()
            << ", predicted impulse force: " << tempEe.estimatedAverageImpulsiveForce.transpose()
            << ", predicted impulse: " << tempEe.estimatedImpulse.transpose() << std::endl;
}
void mi_qpEstimator::print() const
{
  /*
  std::cout<<"The QP estimator params are: "<<std::endl<<"Dim: " <<getImpactModel()->getDim() <<", Dof: "<<getDof()<<",
  coeR: "<<getImpactModel()->getCoeRes()<<", coeF: "<<getImpactModel()->getCoeFricDe()<<", impact duration:
  "<<getImpactModel()->getImpactDuration()<<". "<<std::endl;
*/
  std::cout << red<<"The QP estimator: "<< getEstimatorParams().name<<" has an OSD model with the end-effectors: "<<cyan;
  for(auto idx = getOsd()->getEes().begin(); idx != getOsd()->getEes().end(); ++idx)
  {
    std::cout << *idx << " ";
  }
  std::cout << reset<< std::endl;

  std::cout <<red <<"The QP estimator: "<< getEstimatorParams().name <<" has an  OSD model  with established contacts: "<< green;
  for(auto idx = getOsd()->getContactEes().begin(); idx != getOsd()->getContactEes().end(); ++idx)
  {
    std::cout << *idx << " ";
  }
  std::cout << reset <<std::endl;
}
void mi_qpEstimator::calcPerturbedWrench_()
{
//for(auto & idx: endEffectors_)
  for(auto & idx: endEffectors_)
  {
    idx.second.perturbedWrench.force().setZero();
    idx.second.perturbedWrench.couple().setZero();
for(auto & ii:getImpactModels())
    {
      // If this end-effector is applying an impact, we continue without calculating
      if(ii.first == idx.first )
        continue;

      sva::PTransformd X_0_ii = getSimRobot().bodyPosW(ii.first);
      sva::PTransformd X_0_idx = getSimRobot().bodyPosW(idx.first);

      // Impact frame w.r.t. the idx end-effector frame.
      sva::PTransformd X_ii_idx = X_0_idx*X_0_ii.inv();

     // sva::PTransformd X_ii_idx = getSimRobot().bodyPosW(ii.first).inv();

      sva::ForceVecd tempWrench;
      tempWrench.force().setZero();
      tempWrench.force() = getEndeffector(ii.first).estimatedAverageImpulsiveForce;
      tempWrench.couple().setZero();

      idx.second.perturbedWrench +=  X_ii_idx.dualMul(tempWrench);

    }

    //std::cout<<"qpEstimator: The converted wrench of "<<idx.first<<" is: "<< idx.second.perturbedWrench.vector().transpose()<<std::endl;

  }
}

void mi_qpEstimator::readEeJacobiansSolution_(const Eigen::VectorXd & solutionVariables)
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
    double inv_dt = 1.0 / (getImpactModels().begin()->second->getImpactDuration());

    idx->second.estimatedImpulse = solutionVariables.segment(getDof() + location, getEstimatorParams().dim);

    idx->second.estimatedAverageImpulsiveForce = idx->second.estimatedImpulse * inv_dt;
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
      idx->second.checkForce.setZero();

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
    } // end of Lagrange Multipliers

    // tempEe.jacobianDeltaF = tempA_dagger_ee;

    // std::cout<<"test "<<std::endl;
    // std::cout<<"test "<<std::endl;

    // std::cout<<"tempJ_T size is: "<<tempJ_T.rows() <<", "<<tempJ_T.cols()<<std::endl;

    // std::cout<<"test "<<std::endl;

    // std::cout<<"test "<<std::endl;

  } // end of the for end-effector loop.

}

} // namespace mc_impact
