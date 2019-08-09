# include "mi_qpEstimator.h"


mi_qpEstimator::mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> & osdPtr,
		const struct qpEstimatorParameter params
		): simRobot_(simRobot), osdPtr_(osdPtr), params_(params)
{
  impactModelPtr_.reset(new mi_impactModel(getSimRobot(), getOsd_(), params_.impactBodyName, params_.impactDuration, params_.coeFrictionDeduction, params_.coeRes, params_.dim));

  eqConstraints_.push_back(std::make_shared<mi_iniEquality>(getOsd_(), getImpactModel().get(), false));

  if(params_.useJsd)
    eqConstraints_.push_back(std::make_shared<mi_jsdEquality>(getOsd_(), getImpactModel()->getImpactBody())); 

  if(params_.useOsd)
    eqConstraints_.push_back(std::make_shared<mi_invOsdEquality>(getOsd_(), getImpactModel().get()));

  std::cout<<"Created QP estimator constraint. "<<std::endl;

  initializeQP_(); 
  std::cout<<"the QP-based impulse estimator is created. "<<std::endl;
}
void mi_qpEstimator::initializeQP_()
{
  numVar_  = getDof() + getImpactModel()->getDim()*getOsd_()->getEeNum();

  numEq_ = 0;
  for (auto idx = eqConstraints_.begin(); idx != eqConstraints_.end(); ++idx)
  {
   numEq_+= (*idx)->nrEq();
  }

  solver_.resize(getNumVar_(), getNumEq_(), Eigen::lssol::QP2);

  xl_.resize(getNumVar_()); xu_.resize(getNumVar_());
  xl_ = xl_.setOnes()*-std::numeric_limits<double>::infinity();
  xu_ = xu_.setOnes()*std::numeric_limits<double>::infinity();

  p_.resize(getNumVar_());
  p_.setZero();

  Q_.resize(getNumVar_(), getNumVar_());
  Q_ = Q_.setIdentity()*getQweight();


  C_.resize(getNumEq_(), getNumVar_());
  C_.setZero();
  cu_.resize(getNumEq_());
  cu_.setZero();
  cl_.resize(getNumEq_());
  cl_.setZero();


  jacobianDeltaAlpha_.resize(getDof(), getDof());
  jacobianDeltaAlpha_.setZero();
  jacobianDeltaTau_.resize(getDof(), getDof());
  jacobianDeltaTau_.setZero();

  A_dagger_.resize(
		  getDof() + getEstimatorParams().dim*(static_cast<int>(getOsd_()->getContactEes().size()) + 1), 
		  3);
  A_dagger_.setZero();

  std::cout<<"Reset LSSOL QP estimator variables. "<<std::endl;
}

void  mi_qpEstimator::solveEqQp_(const Eigen::MatrixXd & Q_,const Eigen::VectorXd & p_, const Eigen::MatrixXd & C_, const Eigen::VectorXd & cu_,  Eigen::VectorXd &solution)
{
 Eigen::MatrixXd kkt;
 int kktDim =  getNumVar_() + static_cast<int>(C_.rows());

 kkt.resize(kktDim, kktDim); 
 kkt.setZero();

 kkt.block(0, 0, getNumVar_(), getNumVar_()).setIdentity(getNumVar_(),getNumVar_());
 //kkt.block(0, 0, getNumVar_(), getNumVar_())=  kkt.block(0, 0, getNumVar_(), getNumVar_())*2;


 kkt.block(getNumVar_(), 0, C_.rows(), C_.cols()) = C_;
 kkt.block(0, getNumVar_(), C_.cols(), C_.rows()) = C_.transpose();

 Eigen::VectorXd b; 
 b.resize(kktDim); b.setZero();
 b.segment(getNumVar_(), cu_.rows() ) = cu_;

 // solve the least squares problem: 
 //Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp_kkt(kkt);
 //solution.resize(kktDim);
 //solution = lu_decomp_kkt.inverse()*cu_;
 auto tempInverse =  kkt.completeOrthogonalDecomposition().pseudoInverse();
 solution = tempInverse*b;

 A_dagger_ = tempInverse.block(0, A_dagger_.cols() - 3, A_dagger_.rows(), A_dagger_.cols());
 /*
 Eigen::MatrixXd kktInverse;
 pseudoInverse_(kkt, kktInverse);
 solution = kktInverse*b; 
 */
}

void mi_qpEstimator::pseudoInverse_(const Eigen::MatrixXd & input, Eigen::MatrixXd & output, double tolerance)
{
  auto svd = input.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);

  const auto &singularValues = svd.singularValues();
  Eigen::MatrixXd singularValuesInv(input.cols(), input.rows());

  singularValuesInv.setZero();

  for (unsigned int i = 0; i < singularValues.size(); ++i) {
	  if (singularValues(i) > tolerance)
	  {
		  singularValuesInv(i, i) = 1.0 / singularValues(i);
	  }
	  else
	  {
		  singularValuesInv(i, i) = 0.0;
	  }
  }

  output = svd.matrixV() * singularValuesInv * svd.matrixU().adjoint();
}
void mi_qpEstimator::update(const Eigen::Vector3d & surfaceNormal)
{
  impactModelPtr_->update(surfaceNormal);
  int count = 0;
  // Update the constraints
  for(auto idx = eqConstraints_.begin(); idx != eqConstraints_.end(); ++idx)
  {
    (*idx)->update();

    C_.block(count , 0, (*idx)->nrEq(), getNumVar_()) = (*idx)->AEq();
    cl_.segment(count, (*idx)->nrEq())  = (*idx)->bEq();
    cu_.segment(count, (*idx)->nrEq())  = (*idx)->bEq();

    count +=  (*idx)->nrEq();
  }

  Eigen::VectorXd solutionVariables; 
  solutionVariables.resize(getNumVar_());
  solutionVariables.setZero();

  
  if(params_.useLagrangeMultiplier){

    Eigen::VectorXd b_temp = cu_.segment(getNumVar_(), C_.rows());
    Eigen::MatrixXd A_temp = C_;
    solveEqQp_(Q_, p_, C_, cu_, solutionVariables); 
    //solutionVariables = A_temp.completeOrthogonalDecomposition().pseudoInverse()*b_temp;

  }else{
    solver_.solve(xl_, xu_, Q_, p_, C_, cl_, cu_ );
    solutionVariables = solver_.result();
  }
  
  
  auto tempJac = getImpactModel()->getProjector()*getOsd_()->getJacobian(getImpactModel()->getImpactBody());
  jointVelJump_ = solutionVariables.segment(0, getDof());

  jacobianDeltaTau_.resize(getDof(),3);
  jacobianDeltaTau_.setZero();
  for(auto idx = getOsd_()->getEes().begin(); idx!=getOsd_()->getEes().end(); ++idx)
  {
    int eeIndex = getOsd_()->nameToIndex_(*idx);
    int location =  getOsd_()->getJacobianDim()*eeIndex;
    auto tempEe =  getEndeffector_(*idx);

    tempEe.estimatedImpulse =  solutionVariables.segment(getDof() + location, getOsd_()->getJacobianDim());
	    
    tempEe.estimatedAverageImpulsiveForce= 
	   tempEe.estimatedImpulse/getImpactModel()->getImpactDuration();

    tempEe.eeVJump = getOsd_()->getJacobian(*idx)*jointVelJump_; 


    auto tempA_dagger_ee = A_dagger_.block(getDof() + location, 0 ,3 ,3);

    tempEe.jacobianDeltaF = tempA_dagger_ee*tempJac;

    jacobianDeltaTau_ +=  getOsd_()->getJacobian(*idx).transpose()*tempA_dagger_ee;
  }

  jacobianDeltaAlpha_ = A_dagger_.block(0, 0, getDof(), 3)*tempJac;
  jacobianDeltaTau_ = jacobianDeltaTau_*tempJac;
}



const endEffector & mi_qpEstimator::getEndeffector( const std::string& name) 
{
   return getEndeffector_(name);
}

endEffector & mi_qpEstimator::getEndeffector_( const std::string& name) 
{
  const auto opt = endEffectors_.find(name);
  if(opt != (endEffectors_.end()))
  {
    return opt->second;
  }
  else
  {
    throw std::runtime_error(std::string("getEndeffector: '") +name 
                             + std::string("' is not found."));
  }

}

void mi_qpEstimator::addEndeffector(std::string eeName)
{
  //addOptVariables_(eeName, getDim());
  Eigen::VectorXd tempForce;
  tempForce.resize(getImpactModel()->getDim());
  tempForce.setZero();

  Eigen::VectorXd tempJ;
  tempJ.resize(getEstimatorParams().dim, getDof());
  tempJ.setZero();

  //optVariables_[name] = {dim, optVariables_.size() };
  endEffectors_[eeName] = {tempForce, tempForce, tempForce, tempJ} ;
  //endEffectorNames_.push_back(eeName);

  if(!getOsd_()->addEndeffector(eeName))
  {
    throw std::runtime_error(std::string("OSD failed to add endeffector! ") + eeName);
  }

}

const Eigen::VectorXd  & mi_qpEstimator::getPredictedImpulse(const std::string & bodyName)  
{
  return getEndeffector(bodyName).estimatedImpulse;
}

void mi_qpEstimator::print() const
{
  std::cout<<"The QP estimator params are: "<<std::endl<<"Dim: " <<getImpactModel()->getDim() <<", Dof: "<<getDof()<<", coeR: "<<getImpactModel()->getCoeRes()<<", coeF: "<<getImpactModel()->getCoeFricDe()<<", impact duration: "<<getImpactModel()->getImpactDuration()<<". "<<std::endl;

 std::cout<<"The QP estimator has end-effectors: "; 
 for (auto idx = getOsd_()->getEes().begin(); idx!=getOsd_()->getEes().end(); ++idx)
 {
   std::cout<<*idx<<" "; 
 }
 std::cout<<std::endl;

 std::cout<<"The QP estimator has end-effectors with established contact : ";
 for (auto idx = getOsd_()->getContactEes().begin(); idx!=getOsd_()->getContactEes().end(); ++idx)
 {
   std::cout<<*idx<<" "; 
 }
 std::cout<<std::endl;

}
