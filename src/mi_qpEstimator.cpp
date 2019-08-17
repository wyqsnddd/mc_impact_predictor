# include "mi_qpEstimator.h"


mi_qpEstimator::mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> osdPtr,
		const struct qpEstimatorParameter params
		): simRobot_(simRobot), osdPtr_(osdPtr), params_(params)
{
  
  if(params_.useJsd)
    eqConstraints_.push_back(std::make_shared<mi_jsdEquality>(getOsd(), getEstimatorParams().impactBodyNames )); 

  if(params_.useOsd)
    eqConstraints_.push_back(std::make_shared<mi_invOsdEquality>(getOsd()));


  for (std::vector<std::string>::const_iterator idx = params.impactBodyNames.begin(); idx!=params.impactBodyNames.end(); ++idx)
  {
   impactModels_[*idx] = std::make_shared<mi_impactModel>(getSimRobot(), getOsd(), *idx, params_.impactDuration, params_.timeStep, params_.coeFrictionDeduction, params_.coeRes, params_.dim);

  //eqConstraints_.push_back(std::make_shared<mi_iniEquality>(getOsd(), getImpactModel(const_cast<std::string&>(*idx)).get(), false));
  eqConstraints_.push_back(std::make_shared<mi_iniEquality>(getOsd(), getImpactModel(*idx).get(), false));
  }
  vector_A_dagger_.resize(params.impactBodyNames.size());

  std::cout<<"Created QP estimator constraint. "<<std::endl;

  initializeQP_(); 
  std::cout<<"the QP-based impulse estimator is created. "<<std::endl;
}

void mi_qpEstimator::initializeQP_()
{
  numVar_  = getDof() + 3*getOsd()->getEeNum();

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

  jointVelJump_.resize(getDof());
  tauJump_.resize(getDof());
  jointVelJump_.setZero();
  tauJump_.setZero();

  jacobianDeltaAlpha_.resize(getDof(), getDof());
  jacobianDeltaAlpha_.setZero();
  jacobianDeltaTau_.resize(getDof(), getDof());
  jacobianDeltaTau_.setZero();

  for (size_t idx = 0; idx< vector_A_dagger_.size(); ++idx)
  {
    vector_A_dagger_[idx].resize(
		    
		  getDof() + getEstimatorParams().dim*(getOsd()->getEeNum()), 3
		    );
    vector_A_dagger_[idx].setZero();
  }
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
 

 //Eigen::MatrixXd tempInverse =  kkt.completeOrthogonalDecomposition().pseudoInverse();

 //Suppose kkt has full column rank, we perform cholesky decomposition to compute pinv = (A^*A)^{-1}A^*
 Eigen::MatrixXd tempInverse =  (kkt.transpose()*kkt).inverse()*kkt.transpose();


 solution = tempInverse*b;
 //solution = kkt.completeOrthogonalDecomposition().pseudoInverse()*b;

 //std::cout<<"tempInverse size is: "<<tempInverse.rows() <<", "<<tempInverse.cols()<<std::endl;
 
  //std::cout<<"test 2" <<std::endl;
 for (size_t idx =0; idx< impactModels_.size(); ++idx)
 {
   vector_A_dagger_[idx] = tempInverse.block(0, tempInverse.cols() - 3*static_cast<int>(idx + 1), vector_A_dagger_[idx].rows(), 3);
 
 }

  //std::cout<<"test 3" <<std::endl;
 //A_dagger_ =  //tempInv_= tempInverse.block(0, tempInverse.cols() - 3, tempInverse.rows(), 3);
 tempInv_= tempInverse;
}

void mi_qpEstimator::updateImpactModels_(const std::map<std::string, Eigen::Vector3d> & surfaceNormals)
{
  for (std::map<std::string, Eigen::Vector3d>::const_iterator idx = surfaceNormals.begin(); idx!= surfaceNormals.end(); ++idx)
  {
    //(getImpactModel(const_cast<std::string & >(idx->first)))->update(idx->second);
    (getImpactModel(idx->first))->update(idx->second);
  }
}
void mi_qpEstimator::update(const std::map<std::string, Eigen::Vector3d> & surfaceNormals)
{
  if(surfaceNormals.size()!=impactModels_.size() ){
    throw std::runtime_error(
		    std::string("mi_qpEstimator-update: surfaceNormals size(") 
		    + std::to_string(static_cast<int>(surfaceNormals.size()))
		    + std::string(") does not match impact predictor impact number (") 
		    + std::to_string(impactModels_.size())
		    + std::string(").")
		    );
  }

  updateImpactModels_(surfaceNormals);

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

  //std::cout<<"test 0" <<std::endl;
  if(params_.useLagrangeMultiplier){

    //Eigen::VectorXd b_temp = cu_.segment(getNumVar_(), C_.rows());
    //Eigen::MatrixXd A_temp = C_;
    solveEqQp_(Q_, p_, C_, cu_, solutionVariables); 
    //solutionVariables = A_temp.completeOrthogonalDecomposition().pseudoInverse()*b_temp;

  }else{
    solver_.solve(xl_, xu_, Q_, p_, C_, cl_, cu_ );
    solutionVariables = solver_.result();
  }
  
  
  //std::cout<<"test 1" <<std::endl;
  //Eigen::MatrixXd tempJac = getImpactModel()->getProjector();
  jointVelJump_ = solutionVariables.segment(0, getDof());
  tauJump_.setZero();

  jacobianDeltaTau_.resize(getDof(), getDof());

  //std::cout<<"jacobianDeltaTau_ size is: "<<jacobianDeltaTau_.rows() <<", "<<jacobianDeltaTau_.cols()<<std::endl;
  jacobianDeltaTau_.setZero();

  

  for(auto idx = getOsd()->getEes().begin(); idx!=getOsd()->getEes().end(); ++idx)
  {
    int eeIndex = getOsd()->nameToIndex_(*idx);
    int location =  getOsd()->getJacobianDim()*eeIndex;
    auto & tempEe =  getEndeffector_(*idx);
    //double inv_dt = 1.0/getImpactModel(const_cast<std::string&>(getOsd()->getContactEes()[0]))->getImpactDuration();
    double inv_dt = 1.0/(getImpactModels().begin()->second->getImpactDuration());

    tempEe.estimatedImpulse =  solutionVariables.segment(getDof() + location, getOsd()->getJacobianDim());
    
    tempEe.estimatedAverageImpulsiveForce= 
	   tempEe.estimatedImpulse*inv_dt;

    tempEe.eeVJump = getOsd()->getJacobian(*idx)*jointVelJump_; 


    //std::cout<<"A_dagger_ size is: "<<A_dagger_.rows() <<", "<<A_dagger_.cols()<<std::endl;
    //std::cout<<"Dof is: "<<getDof()<<std::endl;
    //std::cout<<"location is: "<<location<<std::endl;


    //std::cout<<"tempA_dagger_ee is: "<<std::endl<<tempA_dagger_ee<<std::endl;
    if(getEstimatorParams().useLagrangeMultiplier){

      //int kktDim =  getNumVar_() + static_cast<int>(C_.rows());
      //Eigen::VectorXd b; 
      //b.resize(kktDim); b.setZero();
      //b.segment(getNumVar_(), cu_.rows() ) = cu_;

      //auto solution = tempInv_*b; 

      //double inv_dt = 1.0/getImpactModel()->getImpactDuration();
      //Eigen::MatrixXd test_A = tempInv_.block(getDof() + location, 0, 3, 3)  ;
      //tempEe.checkForce = test_A*getImpactModel()->getEeVelocityJump()*inv_dt;
      //tempEe.checkForce = inv_dt*solution.segment(getDof() + location, 3);

      tempEe.jacobianDeltaF.resize(3, getDof());
      tempEe.jacobianDeltaF.setZero();
      tempEe.checkForce.setZero();

      //for (int i?i = 0; ii< static_cast<int>(impactModels_.size()); ++ii )
      unsigned int iiA = 0;
      Eigen::MatrixXd tempJ_T = getOsd()->getJacobian(*idx).transpose();
      for (auto impactIdx = impactModels_.begin(); impactIdx != impactModels_.end(); ++impactIdx, ++iiA)
      {
      
        Eigen::Matrix3d tempA_dagger_ee = vector_A_dagger_[iiA].block(getDof() + location, 0 ,3 ,3);
      
        tempEe.jacobianDeltaF += tempA_dagger_ee*impactIdx->second->getProjector();
	
        tempEe.checkForce += inv_dt*tempEe.jacobianDeltaF*impactIdx->second->getJointVel();

	jacobianDeltaTau_ +=  tempJ_T*tempA_dagger_ee*impactIdx->second->getProjector();
      }

      tauJump_ += tempJ_T*tempEe.estimatedAverageImpulsiveForce; 


      //tempEe.checkForce = inv_dt* tempInv_.block(getDof() + location, tempInv_.cols() - 3, 3, 3)*getImpactModel()->getEeVelocityJump();

      //tempEe.checkForce = inv_dt* tempInv_.block(getDof() + location, tempInv_.cols() , 3, tempInv_.cols())*b;
      //tempEe.checkForce(0) = (b.segment(0, b.rows()-3)).norm();
      //std::cout<<"b is: " <<std::endl<<b.transpose()<<std::endl;
/*
      tempEe.checkForce(0) = (
		      //jointVelJump_  -  tempInv_.block(0, tempInv_.cols() - 3, getDof(), 3)* getImpactModel()->getEeVelocityJump()
		      //jointVelJump_  -  tempInv_.block(0, tempInv_.cols() - 3, getDof(), 3)* tempJac* getImpactModel()->getJointVel()
		      jointVelJump_  -  A_dagger_.block(0, 0, getDof(), 3)* tempJac* getImpactModel()->getJointVel()
		      ).norm();
*/
      
      //std::cout<<"the impulse difference is: "<<(tempEe.estimatedImpulse -  tempA_dagger_ee*getImpactModel()->getEeVelocityJump()).norm()<<std::endl;
    } // end of Lagrange Multipliers

    //tempEe.jacobianDeltaF = tempA_dagger_ee;

    //std::cout<<"test "<<std::endl;
    //std::cout<<"test "<<std::endl;

    //std::cout<<"tempJ_T size is: "<<tempJ_T.rows() <<", "<<tempJ_T.cols()<<std::endl;
    
        //std::cout<<"test "<<std::endl;

    //std::cout<<"test "<<std::endl;
   
  }// end of the for end-effector loop. 

   
  if(getEstimatorParams().useLagrangeMultiplier){
  unsigned int iiA = 0;
  for (auto impactIdx = impactModels_.begin(); impactIdx != impactModels_.end(); ++impactIdx, ++iiA)
  {
    jacobianDeltaAlpha_ = vector_A_dagger_[iiA].block(0, 0, getDof(), 3)*impactIdx->second->getProjector();
  }
  }

/*
  if(getEstimatorParams().useLagrangeMultiplier)
    //jacobianDeltaAlpha_ = A_dagger_.block(0, 0, getDof(), 3)*tempJac;
    jacobianDeltaAlpha_ = tempInv_.block(0, tempInv_.cols() - 3, getDof(), 3)*tempJac;
*/
    //jacobianDeltaAlpha_ = A_dagger_.block(0, 0, getDof(), 3);
  //jacobianDeltaTau_ = jacobianDeltaTau_*tempJac;
}



const endEffector & mi_qpEstimator::getEndeffector( const std::string& name) 
{
   return getEndeffector_(name);
}

endEffector & mi_qpEstimator::getEndeffector_( const std::string& name) 
{
  auto opt = endEffectors_.find(name);
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

const std::shared_ptr<mi_impactModel> & mi_qpEstimator::getImpactModel(const std::string & eeName) 
{
  auto idx = impactModels_.find(eeName);
  if(idx!= (impactModels_.end()))
  {
    return idx->second;
  }
  else
  {
    throw std::runtime_error(std::string("getImpactModel: '") +eeName 
                             + std::string("' is not found."));
  }
}
void mi_qpEstimator::addEndeffector(std::string eeName)
{
  //addOptVariables_(eeName, getDim());
  Eigen::Vector3d tempForce;
  //tempForce.resize(getImpactModel()->getDim());
  tempForce.setZero();

  Eigen::VectorXd tempJ;
  tempJ.resize(getEstimatorParams().dim, getDof());
  tempJ.setZero();

  //optVariables_[name] = {dim, optVariables_.size() };
  endEffectors_[eeName] = {tempForce, tempForce, tempForce, tempForce, tempJ} ;
  //endEffectorNames_.push_back(eeName);

  if(!getOsd()->addEndeffector(eeName))
  {
    throw std::runtime_error(std::string("OSD failed to add endeffector! ") + eeName);
  }

}

void mi_qpEstimator::print(const std::string & eeName)
{
  auto & tempEe = getEndeffector(eeName);
  std::cout<<"Endeffector "<<eeName<<" predicted vel jump: "<<tempEe.eeVJump.transpose()<<", predicted impulse force: "<<tempEe.estimatedAverageImpulsiveForce.transpose()<<", predicted impulse: "<<tempEe.estimatedImpulse.transpose()<<std::endl;
}
void mi_qpEstimator::print() const
{
	/*
  std::cout<<"The QP estimator params are: "<<std::endl<<"Dim: " <<getImpactModel()->getDim() <<", Dof: "<<getDof()<<", coeR: "<<getImpactModel()->getCoeRes()<<", coeF: "<<getImpactModel()->getCoeFricDe()<<", impact duration: "<<getImpactModel()->getImpactDuration()<<". "<<std::endl;
*/
 std::cout<<"The QP estimator has end-effectors: "; 
 for (auto idx = getOsd()->getEes().begin(); idx!=getOsd()->getEes().end(); ++idx)
 {
   std::cout<<*idx<<" "; 
 }
 std::cout<<std::endl;

 std::cout<<"The QP estimator has end-effectors with established contact : ";
 for (auto idx = getOsd()->getContactEes().begin(); idx!=getOsd()->getContactEes().end(); ++idx)
 {
   std::cout<<*idx<<" "; 
 }
 std::cout<<std::endl;

}
