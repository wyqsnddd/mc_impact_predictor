# include "mi_qpEstimator.h"


mi_qpEstimator::mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> & osdPtr,
		const struct qpEstimatorParameter params
		): simRobot_(simRobot), osdPtr_(osdPtr), params_(params)
{
  impactModelPtr_.reset(new mi_impactModel(getSimRobot(), getOsd_(), params_.impactBodyName, params_.impactDuration, params_.coeFrictionDeduction, params_.coeRes, params_.dim));
  eqConstraints_.push_back(std::make_shared<mi_iniEquality>(getOsd_(), getImpactModel().get()));
  eqConstraints_.push_back(std::make_shared<mi_jsdEquality>(getOsd_()));
  eqConstraints_.push_back(std::make_shared<mi_osdEquality>(getOsd_()));

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
  std::cout<<"Reset LSSOL QP estimator variables. "<<std::endl;
}

void mi_qpEstimator::update(const Eigen::Vector3d & surfaceNormal)
{
  impactModelPtr_->update(surfaceNormal);
  int count = 0;
  for(auto idx = eqConstraints_.begin(); idx != eqConstraints_.end(); ++idx)
  {
    (*idx)->update();

    C_.block(count , 0, (*idx)->nrEq(), getNumVar_()) = (*idx)->AEq();
    cl_.segment(count, (*idx)->nrEq())  = (*idx)->bEq();
    cu_.segment(count, (*idx)->nrEq())  = (*idx)->bEq();

    count +=  (*idx)->nrEq();
  }

  solver_.solve(xl_, xu_, Q_, p_, C_, cl_, cu_ );
  
  Eigen::VectorXd solutionVariables = solver_.result();

  jointVelJump_ = solutionVariables.segment(0, getDof());

  for(auto idx = getOsd_()->getEes().begin(); idx!=getOsd_()->getEes().end(); ++idx)
  {
    int eeIndex = getOsd_()->nameToIndex_(*idx);
    int location =  getOsd_()->getJacobianDim()*eeIndex;
    getEndeffector_(*idx).estimatedImpulse =  solutionVariables.segment(getDof() + location, getOsd_()->getJacobianDim());
	    
    getEndeffector_(*idx).estimatedAverageImpulsiveForce= 
	   getEndeffector(*idx).estimatedImpulse/getImpactModel()->getImpactDuration();
      
  }
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

  //optVariables_[name] = {dim, optVariables_.size() };
  endEffectors_[eeName] = {endEffectors_.size(), tempForce, tempForce} ;
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

void mi_qpEstimator::print()
{
  std::cout<<"The QP estimator params are: "<<std::endl<<"Dim: " <<getImpactModel()->getDim() <<", Dof: "<<getDof()<<", coeR: "<<getImpactModel()->getCoeRes()<<", coeF: "<<getImpactModel()->getCoeFricDe()<<", impact duration: "<<getImpactModel()->getImpactDuration()<<". "<<std::endl;

}
