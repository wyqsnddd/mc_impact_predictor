# include "mi_qpEstimator.h"


mi_qpEstimator::mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> & osdPtr,
		const struct qpEstimatorParameter params
		//const std::string & impactBodyName,
		//int  dim,
	  	//const std::string & solverName,
		//double convergenceThreshold
		): simRobot_(simRobot), osdPtr_(osdPtr), params_(params)
	//impactBodyName_(impactBodyName), dim_(dim), solverName_(solverName), convergenceThreshold_(convergenceThreshold)
{
  jointVelJump_.resize(getOsd_()->getDof());
  jointVelJump_.setZero();
  // addOptVariables_(std::string("deltaQDot"), osdPtr_->getDof());
  std::cout<<"the QP-based impulse estimator is created. "<<std::endl;
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
   
  std::vector<double> solutionVariables = solver_.solveQP(this);

  // Fill the end-effector data structure 
  
  jointVelJump_ = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(&solutionVariables[0], getOsd_()->getDof());

  int contactCounter = 0; 
  std::vector<std::string> cEes = getOsd_()->getContactEes();
  for ( auto idx = cEes.begin(); idx != cEes.end(); ++idx)
  {
   contactCounter++;
   Eigen::VectorXd temp = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(&solutionVariables[contactCounter*getDim()], getDim());
   getEndeffector_(*idx).estimatedImpulse = temp;
   getEndeffector_(*idx).estimatedAverageImpulsiveForce= temp/getImpactDuration_();
   //predictedContactForce_[*idx] = temp; 
 }

}



/*
void mi_qpEstimator::addOptVariables_( const std::string& name, const unsigned & dim)
{
  optVariables_[name] = {dim, optVariables_.size() };
}
*/

endEffector & mi_qpEstimator::getEndeffector_( const std::string& name) 
{
  std::map<std::string, endEffector>::iterator opt = endEffectors_.find(name);

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


const endEffector & mi_qpEstimator::getEndeffector( const std::string& name) const
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
  endEffectors_[eeName] = {endEffectors_.size(), tempForce, tempForce} ;

  if(!getOsd_()->addEndeffector(eeName))
  {
    throw std::runtime_error(std::string("OSD failed to add endeffector! ") + eeName);
  }

}

const Eigen::VectorXd  & mi_qpEstimator::getPredictedImpulse(const std::string & bodyName) const
{
  return getEndeffector(bodyName).estimatedImpulse;
}


void qp_solver::jsdImpulseConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data)
{

  const mi_qpEstimator * estimatorPtr = static_cast<mi_qpEstimator* >(f_data);
  // Preparation: 
  //constraintDim = static_cast<unsigned> (estimatorPtr->getDof());
  //stateDim = static_cast<k>
  //result = new double [];


  //const unsigned dim(constraintDataPtr->estimatorPtr->getOptVariablesInfo_("deltaQDot").dim);
  //size_t startingNumber = constraintDataPtr->estimatorPtr->getOptVariablesInfo_("deltaQDot").startingNumber;

  Eigen::Map<const Eigen::VectorXd> state(x, stateDim);
  int dim = estimatorPtr->getDim();
  //std::cout<<"jsd: the dim is: "<<dim<<std::endl;
  int dof = estimatorPtr->getDof();

  Eigen::VectorXd temp;
  temp.resize(dim);

  for(auto idx = estimatorPtr->endEffectors_.begin(); idx!=estimatorPtr->endEffectors_.end(); idx++)
  {
    temp += estimatorPtr->getOsd_()->getJacobian(idx->first)*state.segment(dof + dim*idx->second.startingNumber, dim);
  }

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


void qp_solver::testEq( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data)
{
 
  const mi_qpEstimator* estimatorPtr =  static_cast<mi_qpEstimator* >(f_data);
  result = new double [2];
  *result = x[0] - 1.0;
  *(result+1) = x[1];
  //*(result+2) = x[2];
  grad = new double [2];
  grad[0] = 1; 
  grad[1] = 1; 
  //grad[2] = 1; 
}
void qp_solver::iniConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data)
{
  const mi_qpEstimator* estimatorPtr =  static_cast<mi_qpEstimator* >(f_data);

  Eigen::Map<const Eigen::VectorXd> state(x, stateDim);
  //int dim = estimatorPtr->getDim();
  int dof = estimatorPtr->getDof();
  //std::cout<<"ConstraintDim is: "<< constraintDim <<std::endl; 
  // Evaluate the result: 
  result = new double [constraintDim];
  Eigen::Map<Eigen::VectorXd> value(result, constraintDim);
  //value = estimatorPtr->getOsd_()->getJacobian(estimatorPtr->getImpactBody())*state.segment(0, dof) - estimatorPtr->getEeVelocityJump();
  //value = state.segment(0, 3) - estimatorPtr->getEeVelocityJump(); 
  //value = state.segment(0, 3) - Eigen::Vector3d::Ones(); 

  value.setZero();
  //value = state.segment(0, 3) - Eigen::Vector3d::Ones(); 
  grad = new double [constraintDim];
  Eigen::Map<Eigen::VectorXd> valueGrad(grad, constraintDim);
  valueGrad.setZero();

}



double qp_solver::objFunction(const std::vector<double> &x, std::vector<double> &grad, void *obj_data)
{
  qpObjData* objDataPtr = static_cast<qpObjData*> (obj_data);

  Eigen::Map<const Eigen::VectorXd> q( x.data(), x.size());

  Eigen::VectorXd gradient = objDataPtr->H*q;

  for (std::size_t  ii = 0; ii <grad.size(); ++ii)
  {
    grad[ii] = gradient(static_cast<int>(ii));
  }


  return (double)(0.5*q.transpose()*objDataPtr->H*q); 

}

std::vector<double> & qp_solver::solveQP(mi_qpEstimator* estimatorPtr)
{
  // (0) Init
  //int numVar = estimatorPtr->getDof() + static_cast<int>(estimatorPtr->getDim()*estimatorPtr->endEffectors_.size() );
  //int numVar =  static_cast<int>(estimatorPtr->getDim()*estimatorPtr->endEffectors_.size() );
  int numVar = 3; 
  std::cout<<"QP:: The number of variables is: "<<numVar<<std::endl;

  //nlopt::opt opt(nlopt::LD_SLSQP, static_cast<unsigned int>(numVar));
  // Round-off: 
  //nlopt::opt opt(nlopt::LD_SLSQP, static_cast<unsigned int>(numVar));
  
  nlopt::opt opt(nlopt::LD_CCSAQ, numVar);
  // Never returns: 
  //nlopt::opt opt(nlopt::LN_COBYLA, static_cast<unsigned int>(numVar));
  //nlopt::opt opt(nlopt::LD_AUGLAG, static_cast<unsigned int>(numVar));
  //nlopt::opt opt(nlopt::LD_AUGLAG_EQ, static_cast<unsigned int>(numVar));
  
  // invalid argument
  // nlopt::opt opt(nlopt::GN_ISRES, static_cast<unsigned int>(numVar));
  //nlopt::opt opt(nlopt::LN_SBPLX, static_cast<unsigned int>(numVar));
  // nlopt::opt opt(nlopt::LD_LBFGS, static_cast<unsigned int>(numVar));
  // nlopt::opt opt(nlopt::LD_MMA, static_cast<unsigned int>(numVar));
  
  // (1) set the objective: 
  qpObjData objData = {Eigen::MatrixXd::Identity(numVar, numVar) };
  void * objDataPtr  = &objData;
  opt.set_min_objective(qp_solver::objFunction, objDataPtr);
  opt.set_xtol_rel(estimatorPtr->getThreshold_());
  // (2) Add the constraints:  
  
  void * conDataPtr  = static_cast<void*>(estimatorPtr);
  std::vector<double> jsdTolerance, osdTolerance, iniTolerance, testTolerance;
  testTolerance.resize(2);
  testTolerance.assign(2, 0.01);
  opt.add_inequality_mconstraint(qp_solver::testEq, conDataPtr, testTolerance);
  

  jsdTolerance.resize(static_cast<size_t>(estimatorPtr->getDof()));
  jsdTolerance.assign(static_cast<size_t>(estimatorPtr->getDof()), estimatorPtr->getThreshold_());
  //std::cout<<"Tolerance size is: "<<jsdTolerance.size()<<", it stated with: "<<jsdTolerance[0]<<std::endl;

  //std::cout<<"about to add the jsd equality constraint"<<std::endl;
  //opt.add_equality_mconstraint(qp_solver::jsdImpulseConstraintFunction, conDataPtr, jsdTolerance);
  //std::cout<<"added the jsd equality constraint"<<std::endl;
  //opt.add
  osdTolerance.resize(static_cast<size_t>(estimatorPtr->getDim())*estimatorPtr->endEffectors_.size());
  osdTolerance.assign(static_cast<size_t>(estimatorPtr->getDim())*estimatorPtr->endEffectors_.size(), estimatorPtr->getThreshold_());
  //opt.add_equality_mconstraint(qp_solver::osdImpulseConstraintFunction, conDataPtr, osdTolerance);
 
  iniTolerance.resize(static_cast<size_t>(estimatorPtr->getDim()));
  iniTolerance.assign(static_cast<size_t>(estimatorPtr->getDim()), 10*estimatorPtr->getThreshold_());
  //opt.add_equality_mconstraint(qp_solver::iniConstraintFunction, conDataPtr, iniTolerance);
  // (*)Set the initial values
  solution.resize(static_cast<unsigned long>(numVar));
  std::fill(solution.begin(), solution.end(), 0.0);
  
  double minf;
  // ------ solve: 
  try{
    result = opt.optimize(solution, minf);

    //std::cout<<"LCP::solved"<<std::endl; 
    Eigen::Map<const Eigen::VectorXd> q_opt( solution.data(), solution.size());
    
     std::cout << "found minimum at f(" << q_opt.transpose()<< ") = "
        << minf << std::endl;
  
  }
  catch(std::exception &e) {
    std::cout << "qp_solver::nlopt failed: " << e.what() << std::endl;
  }
  // If the result is not in this range, then it failed
  if ( !(nlopt::SUCCESS <= result && result <= nlopt::XTOL_REACHED) )
   throw std::runtime_error("nlopt failed to solve the problem");

  return solution;
}
void mi_qpEstimator::print()
{
  std::cout<<"The QP estimator params are: "<<std::endl<<"Dim: " <<getDim() <<", Dof: "<<getDof()<<", coeR: "<<getCoeRes_()<<", coeF: "<<getCoeFricDe_()<<", soverName: "<<getSolver_()<<", impact duration: "<<getImpactDuration_()<<", convergence threshold: "<<getThreshold_()<<std::endl;

}
