# include "mi_qpEstimator.h"


mi_qpEstimator::mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> & osdPtr,
		const struct qpEstimatorParameter params
		//const std::string & impactBodyName,
		//int  dim,
	  	//const std::string & solverName,
		//double convergenceThreshold
		): simRobot_(simRobot_), osdPtr_(osdPtr), params_(params)
	//impactBodyName_(impactBodyName), dim_(dim), solverName_(solverName), convergenceThreshold_(convergenceThreshold)
{

  // addOptVariables_(std::string("deltaQDot"), osdPtr_->getDof());
  std::cout<<"the QP-based impulse estimatos is created. "<<std::endl;
}

void mi_qpEstimator::update(const Eigen::Vector3d & surfaceNormal)
{
  Eigen::Matrix3d tempProjector = surfaceNormal * surfaceNormal.transpose();
  Eigen::Matrix3d tempNullProjector = Eigen::Matrix3d::Identity() - surfaceNormal * surfaceNormal.transpose();


  Eigen::Matrix3d tempReductionProjector = -((1 + getCoeRes_()) * tempProjector + getCoeFricDe_() * tempNullProjector);
  Eigen::VectorXd alpha = rbd::dofToVector(getSimRobot().mb(), getSimRobot().mbc().alpha);
  Eigen::VectorXd alphaD = rbd::dofToVector(getSimRobot().mb(), getSimRobot().mbc().alphaD);

  deltaV_ = tempReductionProjector * getOsd_()->getJacobian(getImpactBody()) * (alpha + alphaD * getImpactDuration_());

  // run the QP estimator
}



/*
void mi_qpEstimator::addOptVariables_( const std::string& name, const unsigned & dim)
{
  optVariables_[name] = {dim, optVariables_.size() };
}
*/

const endEffector & mi_qpEstimator::getEndeffector_( const std::string& name) const
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
  tempForce.resize(getDim());
  tempForce.setZero();

  //optVariables_[name] = {dim, optVariables_.size() };
  endEffectors_[eeName] = {endEffectors_.size(), tempForce} ;

  if(!getOsd_()->addEndeffector(eeName))
  {
    throw std::runtime_error(std::string("OSD failed to add endeffector! ") + eeName);
  }

}

const Eigen::VectorXd  & mi_qpEstimator::getPredictedImpulse(const std::string & bodyName) const
{
  return getEndeffector_(bodyName).estimatedImpulse;
}


void qp_solver::jsdImpulseConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data)
{
  const mi_qpEstimator * estimatorPtr = static_cast<mi_qpEstimator* >(f_data);

  //const unsigned dim(constraintDataPtr->estimatorPtr->getOptVariablesInfo_("deltaQDot").dim);
  //size_t startingNumber = constraintDataPtr->estimatorPtr->getOptVariablesInfo_("deltaQDot").startingNumber;

  Eigen::Map<const Eigen::VectorXd> state(x, stateDim);
  int dim = estimatorPtr->getDim();
  int dof = estimatorPtr->getDof();

  Eigen::VectorXd temp;
  temp.resize(dim);

  for(auto idx = estimatorPtr->endEffectors_.begin(); idx!=estimatorPtr->endEffectors_.end(); idx++)
  {
    temp += estimatorPtr->getOsd_()->getJacobian(idx->first)*state.segment(dof + dim*idx->second.startingNumber, dim);
  }

  // Evaluate the result: 
  Eigen::Map<Eigen::VectorXd> value(result, constraintDim);
  value = estimatorPtr->getOsd_()->getMassMatrix()*state.segment(0, dof) - temp;
  // Evaluate the gradient: leave the gradient vacant for now. 

}

void qp_solver::osdImpulseConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data)
{
  const mi_qpEstimator* estimatorPtr = static_cast<mi_qpEstimator* >(f_data); 


  Eigen::Map<const Eigen::VectorXd> state(x, stateDim);
  int dim = estimatorPtr->getDim();
  int dof = estimatorPtr->getDof();

  // Get the corresponding impulse 
  //size_t startingNumber = estimatorPtr->getEndeffector_((constraintDataPtr->eeName)).startingNumber;
  // 
  // Evaluate the result: 
  Eigen::Map<Eigen::VectorXd> value(result, constraintDim);
  for(auto idx = estimatorPtr->endEffectors_.begin(); idx!=estimatorPtr->endEffectors_.end(); idx++)
  {

    value.segment(dim*idx->second.startingNumber, dim) = 
	    estimatorPtr->getOsd_()->getEffectiveLambdaMatrix( idx->first)*state.segment(0, dof) - state.segment(dof + dim*idx->second.startingNumber, dim);
	   // *estimatorPtr->getOsd_()->getJacobian(idx->first);
  }

//  value = - state.segment(dof + dim*startingNumber, dim);

}



void qp_solver::iniConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data)
{
  const mi_qpEstimator* estimatorPtr =  static_cast<mi_qpEstimator* >(f_data);

  Eigen::Map<const Eigen::VectorXd> state(x, stateDim);
  //int dim = estimatorPtr->getDim();
  int dof = estimatorPtr->getDof();
   
  // Evaluate the result: 
  Eigen::Map<Eigen::VectorXd> value(result, constraintDim);
  value = estimatorPtr->getOsd_()->getJacobian(estimatorPtr->getImpactBody())*state.segment(0, dof) - estimatorPtr->getEeVelocityJump();

}



double qp_solver::objFunction(const std::vector<double> &x, std::vector<double> &grad, void *obj_data)
{
  quadraticObjData* objDataPtr = static_cast<quadraticObjData*> (obj_data);

  Eigen::Map<const Eigen::VectorXd> q( x.data(), x.size());

  Eigen::VectorXd gradient = objDataPtr->H*q;

  for (std::size_t  ii = 0; ii <grad.size(); ++ii)
  {
    grad[ii] = gradient(static_cast<int>(ii));
  }


  return (double)(0.5*q.transpose()*objDataPtr->H*q); 

}

std::vector<double> & qp_solver::solveQP(const mi_qpEstimator* estimatorPtr)
{
  // (0) Init
  int numVar = static_cast<int>(estimatorPtr->getDof() + estimatorPtr->endEffectors_.size() );
  //std::cout<<"LCP:: The number of variables is: "<<numVar<<std::endl;
  //nlopt::opt opt(estimatorPtr->getSolver_().c_str(), numVar);

  nlopt::opt opt(nlopt::LD_CCSAQ, numVar);

  // (1) set the objective: 
  quadraticObjData objData = {Eigen::MatrixXd::Identity(numVar, numVar) };
  void * objDataPtr  = &objData;
  opt.set_min_objective(qp_solver::objFunction, objDataPtr);
  opt.set_xtol_rel(estimatorPtr->getThreshold_());
  // (2) Add the constraints:  
  std::vector<double> jsdTolerance, osdTolerance, iniTolerance;
  jsdTolerance.resize(estimatorPtr->getDof());
  jsdTolerance.assign(estimatorPtr->getDof(), estimatorPtr->getThreshold_());
  void * jsdDataPtr  = &estimatorPtr;
  opt.add_equality_mconstraint(qp_solver::jsdImpulseConstraintFunction, jsdDataPtr, jsdTolerance);
   
  osdTolerance.resize(estimatorPtr->getDim()*estimatorPtr->endEffectors_.size());
  osdTolerance.assign(estimatorPtr->getDim()*estimatorPtr->endEffectors_.size(), estimatorPtr->getThreshold_());
  void * osdDataPtr = &estimatorPtr; 
  opt.add_equality_mconstraint(qp_solver::osdImpulseConstraintFunction, osdDataPtr, osdTolerance);
 
  iniTolerance.resize(estimatorPtr->getDof());
  iniTolerance.assign(estimatorPtr->getDof(), estimatorPtr->getThreshold_());
  void * iniDataPtr = &estimatorPtr; 
  opt.add_equality_mconstraint(qp_solver::iniConstraintFunction, iniDataPtr, iniTolerance);

  // (*)Set the initial values
  solution.resize(numVar);
  std::fill(solution.begin(), solution.end(), 0.0);
  
  double minf;
  // ------ solve: 
  try{
    result = opt.optimize(solution, minf);

    //std::cout<<"LCP::solved"<<std::endl; 
    Eigen::Map<const Eigen::VectorXd> q_opt( solution.data(), solution.size());
    /*
     std::cout << "found minimum at f(" << q_opt.transpose()<< ") = "
        << minf << std::endl;
  */
  }
  catch(std::exception &e) {
    std::cout << "lcp_solver::nlopt failed: " << e.what() << std::endl;
  }
  // If the result is not in this range, then it failed
  if ( !(nlopt::SUCCESS <= result && result <= nlopt::XTOL_REACHED) )
   throw std::runtime_error("nlopt failed to solve the problem");

  return solution;
}

