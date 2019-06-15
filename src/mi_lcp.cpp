# include "mi_lcp.h"


mi_lcp::mi_lcp(mc_rbdyn::Robot & robot,
		const std::shared_ptr<mi_osd> & osdPtr,
	int  dim)
: robot_(robot), osdPtr_(osdPtr), dim_(dim){
}

void mi_lcp::update_(const Eigen::MatrixXd & Jacobian, const Eigen::MatrixXd & JacobianDot)
{
// (1) Update beta and d 

 Eigen::VectorXd alpha = rbd::dofToVector(getRobot().mb(), getRobot().mbc().alpha);
 Eigen::MatrixXd tempMInv = Jacobian*osdPtr_->getInvMassMatrix();
 Eigen::VectorXd tau = rbd::dofToVector(getRobot().mb(), getRobot().mbc().jointTorque);
 std::cout<<"Tau is: "<<tau.transpose()<<std::endl;
 beta_ = JacobianDot*alpha - tempMInv*osdPtr_->getFD()->C();
 //std::cout<<"C is: "<<osdPtr_->getFD()->C().transpose()<<std::endl;
 d_ = tempMInv*tau + beta_;


 // (2) construct the QP and solve the LCP  
 Eigen::MatrixXd tempLambdaInv = tempMInv*Jacobian.transpose(); 
 
 // (3) solve the lcp problem 
 std::vector<double>testForce = solver_.solveLCP(tempLambdaInv, d_);

 Eigen::Map<const Eigen::VectorXd> q_opt( testForce.data(), testForce.size());
 std::cout<<"LCP says that the contact Force is: "<<q_opt.transpose()<<std::endl;

}

void mi_lcp::update()
{
// We restrict the calculation to the surface normal direction
// The number of contact is: 
 std::size_t contactNum = osdPtr_->getContactNum();
 int dof = osdPtr_->getDof();
 Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(contactNum*getDim_(), dof); 
 Eigen::MatrixXd JacobianDot = Eigen::MatrixXd::Zero(contactNum*getDim_(), dof); 

 int contactCounter = 0; 
 //for ( std::vector<std::string>::iterator idx = osdPtr_->getContactEes().begin(); idx != osdPtr_->getContactEes.end(); ++idx)
 std::vector<std::string> cEes = osdPtr_->getContactEes();
 for ( auto idx = cEes.begin(); idx != cEes.end(); ++idx)
 {
  Jacobian.block(contactCounter*getDim_(), 0, getDim_(), dof) = osdPtr_->getJacobian(*idx);
  JacobianDot.block(contactCounter*getDim_(), 0, getDim_(), dof) = osdPtr_->getJacobianDot(*idx);
  contactCounter++;
 }
 update_(Jacobian, JacobianDot);

}
void mi_lcp::update(std::map<std::string, Eigen::Vector3d> contactSurfaceNormals)
{
// We restrict the calculation to the surface normal direction
// The number of contact is: 
 std::size_t contactNum = osdPtr_->getContactNum();
 int dof = osdPtr_->getDof();
 Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(contactNum, dof); 
 Eigen::MatrixXd JacobianDot = Eigen::MatrixXd::Zero(contactNum, dof); 

 int contactCounter = 0; 
 //for ( std::vector<std::string>::iterator idx = osdPtr_->getContactEes().begin(); idx != osdPtr_->getContactEes.end(); ++idx)
 std::vector<std::string> cEes = osdPtr_->getContactEes();
 for ( auto idx = cEes.begin(); idx != cEes.end(); ++idx)
 {
  auto ee = contactSurfaceNormals.find(*idx);
  
  Jacobian.block(contactCounter*getDim_(), 0, getDim_(), dof) = ee->second.transpose()*osdPtr_->getJacobian(*idx);
  JacobianDot.block(contactCounter*getDim_(), 0, getDim_(), dof) = ee->second.transpose()*osdPtr_->getJacobianDot(*idx);
  contactCounter++;
 }
 update_(Jacobian, JacobianDot);
}



std::vector<double> & lcp_solver::solveLCP(const Eigen::MatrixXd & H, const Eigen::VectorXd & d )
{
  int numVar = static_cast<int>(d.size());
  std::cout<<"LCP:: The number of variables is: "<<numVar<<std::endl;
  nlopt::opt opt(nlopt::LD_CCSAQ, numVar);
  std::vector<double> lb(numVar, 0.0);
  opt.set_lower_bounds(lb);
  quadraticObjData objData = {H, d};
  void * objDataPtr  = &objData;
  opt.set_min_objective(lcp_solver::objFunction, objDataPtr);
  std::cout<<"LCP::objective is set"<<std::endl; 
  opt.set_xtol_rel(1e-2);
  solution.resize(numVar);
  std::fill(solution.begin(), solution.end(), 100.0);
  double minf;
  try{
    result = opt.optimize(solution, minf);

    std::cout<<"LCP::solved"<<std::endl; 
    Eigen::Map<const Eigen::VectorXd> q_opt( solution.data(), solution.size());
    std::cout << "found minimum at f(" << q_opt.transpose()<< ") = "
        << minf << std::endl;
  }
  catch(std::exception &e) {
    std::cout << "lcp_solver::nlopt failed: " << e.what() << std::endl;
  }

  // If the result is not in this range, then it failed
  if ( !(nlopt::SUCCESS <= result && result <= nlopt::XTOL_REACHED) )
   throw std::runtime_error("nlopt failed to solve the problem");
  
  //Eigen::Map<Eigen::VectorXd> solution_eigen( solution.data(), solution.size());
  //std::cout<<"LCP::solution mapped"<<std::endl; 
  //std::cout<<solution_eigen.transpose()<<std::endl;
  return solution;
}

double lcp_solver::objFunction(const std::vector<double> &x, std::vector<double> &grad, void *obj_data)
{
  quadraticObjData* objDataPtr = static_cast<quadraticObjData*> (obj_data);

  std::cout<<"LCP objective:: H is: "<<objDataPtr->H<<std::endl; 
  std::cout<<"LCP objective:: d is: "<<objDataPtr->d.transpose()<<std::endl; 
  Eigen::Map<const Eigen::VectorXd> q( x.data(), objDataPtr->d.size());

  //double temp =  (double)0.5*q.transpose()*objDataPtr->H*q;
  return (double)(0.5*q.transpose()*objDataPtr->H*q) + objDataPtr->d.transpose()*q;

}



