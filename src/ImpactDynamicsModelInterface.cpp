#include "ImpactDynamics/ImpactDynamicsModelInterface.h" 

namespace mc_impact 
{


ImpactDynamicsModel::ImpactDynamicsModel(
		const std::shared_ptr<RobotInterface::Robot> robotPtr,
		const ImpactModelParams & params)
: robotPtr_(robotPtr), params_(params)
{
  std::cout << RoboticsUtils::info << "ImpactDynamicsModel is created." << RoboticsUtils::reset << std::endl;
}

TwoDimModelBridge::TwoDimModelBridge(
		  const std::shared_ptr<RobotInterface::Robot> robotPtr,
                  const ImpactModelParams & params,
		  const TwoDimModelBridgeParams & bridgeParams)
: ImpactDynamicsModel(robotPtr, params), bridgeParams_(bridgeParams)
{

  // Initialize the two-dim-model
  // Energetic coefficient of restitution
  piParams_.e = 0.2;
  // Coefficient of friction
  piParams_.miu = 0.3;
  rotation_.setZero();
  rotationFull_.setIdentity();

  twoDimModelPtr_.reset(new FIDynamics::TwoDimModel(piParams_));

  if(getTwoDimModelBridgeParams().gradientApproximation)
  {
    initializeGradientApproximation_();
  }

  std::cout << RoboticsUtils::info << "TwoDimModelBridge is created." << RoboticsUtils::reset << std::endl;
}


void TwoDimModelBridge::update(const Eigen::Vector3d & impactNormal, const Eigen::Vector3d & impactLinearVel)
{

  // std::cout<<"Updating TwoDimModelBridge"<<std::endl;
  // (0) Compute the whole-body inertia and average velocity
  Eigen::Matrix6d centroidalInertia;
  centroidalInertia.setIdentity();

  // Eigen::Vector6d cm;
  Eigen::Vector6d av;

  getRobot()->computeCentroidalInertia(centroidalInertia, av);

  // std::cout<<"inertia updated"<<std::endl;
  // Assert that the average com velocity is equal to the com velocity
  /*
  assert(mc_impact::areSame(av(3), getRobot().comVelocity()(0)));
  assert(mc_impact::areSame(av(4), getRobot().comVelocity()(1)));
  assert(mc_impact::areSame(av(5), getRobot().comVelocity()(2)));
  */

  rAverageAngularVel_ = av.segment<3>(0);
  rAverageLinearVel_ = av.segment<3>(3);
  rCentroidalInertia_ = centroidalInertia.block<3, 3>(0, 0);

  // std::cout<<"The average angular velocity is: "<<rAverageAngularVel_.transpose()<<std::endl;
  // std::cout<<"The average linear velocity is: "<<getRobot().comVelocity().transpose()<<std::endl;

   // clang-format off
// This is a sample centroidal inertia:
/*
   6.50171    0.0572213     0.282871   1.1147e-17  7.68482e-16 -3.81639e-16
   0.0572213      6.06731       0.1157 -2.58474e-16  3.37024e-17 -3.95517e-16
   0.282871       0.1157      1.81269  1.94289e-16  1.78677e-16 -1.89735e-17
   1.8384e-17  2.13371e-16   3.1225e-17      39.0255  9.96488e-17 -1.82363e-16
   2.75821e-16  5.28301e-17  1.17961e-16  4.28319e-17      39.0255  1.39119e-16
  -5.89806e-16 -4.02456e-16 -1.10453e-18 -6.34367e-16  2.32848e-16      39.0255
*/
  // clang-format on
  // Use the impact body translation
  //sva::PTransformd X_0_ee = getRobot()->bodyPosW(getParams().iBodyName);
  auto type = getRobot()->getImplementationType();

  sva::PTransformd X_0_ee;

  if(type == "McRobot")
  {
     X_0_ee =  getRobot()->bodyPosW(getParams().iBodyName);
  }else if(type == "DartRobot"){
     X_0_ee =  getRobot()->bodyPosW(getParams().iSurfaceName);
  }else
  {
    RoboticsUtils::throw_runtime_error("", __FILE__, __LINE__) ;
  }
  // std::cout<<"vc parameters updated"<<std::endl;

  params_.eePosition = X_0_ee.translation();


  // (3) Update the twoDim model
  updatePiParams_(impactNormal, getParams().eePosition, impactLinearVel);

  twoDimModelPtr_->updateParams(getPlanarImpactParams());
  // std::cout<<"twoDimModelPtr_->updated params "<<std::endl;
  twoDimModelPtr_->update();

  // std::cout<<"twoDimModelPtr_->updated "<<std::endl;
  // (4) Convert the twoDim model solution back to 3D:
  planarSolutionTo3D_();

  /*
  // std::cout<<"converted solution to 3d"<<std::endl;
  if(getTwoDimModelBridgeParams().gradientApproximation)
  {
    // velCases_.clear();

    Eigen::Vector3d contactVel = impactLinearVel;

    double contactVelGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
    double comVelJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];

    int ii = 0;
    // Loop over the caurse grids:
    for(const auto & vel : caurseContactVelocityGrids_)
    {

      contactVel.x() = vel;
      contactVelGrids[ii] = vel;

      if(vel != impactLinearVel.x())
      {
        // Compute if it is not the nominal contact velocity

        updatePiParams_(impactNormal, getParams().eePosition, contactVel);
        // Reset the x component of the impact velocity
        twoDimModelPtr_->updateParams(getPlanarImpactParams());
        // std::cout<<"twoDimModelPtr_->updated params "<<std::endl;
        twoDimModelPtr_->update();
      }

      planarSolutionTo3DPushWall_(velCases_[vel]);

      // velCases_.insert(std::make_pair(vel, getRobotPostImpactStates().impulse));
      // velCases_[vel] = getRobotPostImpactStates().impulse;

      // velCases_[vel] = getRobotPostImpactStates();
      comVelJumpGrids[ii] = velCases_[vel].linearVelJump.x();

      // getRobotPostImpactStates().impulse.x(),  getRobotPostImpactStates().impulse.y(),
      // getRobotPostImpactStates().impulse.z();

      // std::cout<<"For contact vel: "<< vel<<", the impulse is: "<<
      // getRobotPostImpactStates().impulse.transpose()<<std::endl;

      ii++;
    }
    // std::cout<<"The number of velCases are: "<< velCases_.size()<<std::endl;

    double c0, c1, cov00, cov01, cov11, sumsq;
    
    // gsl_fit_linear(contactVelGrids, 1, comVelJumpGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		   
    c0 = 0.0;
    gsl_fit_mul(contactVelGrids, 1, comVelJumpGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &c1, &cov11, &sumsq);

    robotPostImpactStates_.c = c1;
  }
  */

} // end of update

/*
void TwoDimModelBridge::update(const Eigen::Vector3d & impactNormal, StandingStabilityParams & params)
{
  computeGradient_(impactNormal, params);

  auto cvj_normal = getRobotPostImpactStates().linearVelJump.normalized();

  mcZMPAreaPtr_->computeZMPBound(-cvj_normal, params.zmpLowerBound, params.zmpUpperBound);

    computeMaxContactVel_(impactNormal, params.zmpLowerDistance, params.zmpUpperDistance, params);
}
*/
void TwoDimModelBridge::gradientApproximationCalc_(
		  const Eigen::Vector3d & impactNormal, 
		  double & c1)
{
   // params.jumpDirection.x() = cx.c1;
 // params.jumpDirection.y() = cy.c1;
 // params.jumpDirection.z() = cz.c1;

  c1 = robotPostImpactStates_.gradient.coe.norm(); 
  // The sign depends on the angle between the "contact vel" and the "post-impact com velocity jump"
  
  double signTest =  impactNormal.transpose()*robotPostImpactStates_.gradient.coe;

  //std::cout<<"Interface: norm: "<<params.c1<<std::endl;
  c1 = -RoboticsUtils::sgn(signTest) * c1;
  //std::cout<<"c1 is:: "<<params.c1<<std::endl;
  robotPostImpactStates_.c = c1;

}


void TwoDimModelBridge::computeGradient(const Eigen::Vector3d & impactNormal, Eigen::Vector3d & jumpDirection, double & c1)
{
  // (0) Compute the whole-body inertia and average velocity
  Eigen::Matrix6d centroidalInertia;
  centroidalInertia.setIdentity();

  // Eigen::Vector6d cm;
  Eigen::Vector6d av;

  getRobot()->computeCentroidalInertia(centroidalInertia, av);
  rAverageAngularVel_ = av.segment<3>(0);
  rAverageLinearVel_ = av.segment<3>(3);
  rCentroidalInertia_ = centroidalInertia.block<3, 3>(0, 0);

  //sva::PTransformd X_0_ee = getRobot()->bodyPosW(getParams().iBodyName);
  //sva::PTransformd X_0_ee = getRobot()->bodyPosW(getParams().iSurfaceName);
  sva::PTransformd X_0_ee; 

  auto type = getRobot()->getImplementationType();

  if(type == "McRobot")
  {
     X_0_ee =  getRobot()->bodyPosW(getParams().iBodyName);
  }else if(type == "DartRobot"){
     X_0_ee =  getRobot()->bodyPosW(getParams().iSurfaceName);
  }else
  {
    RoboticsUtils::throw_runtime_error("", __FILE__, __LINE__) ;
  }
  params_.eePosition = X_0_ee.translation();

  // (1) Approximate the gradient. 
  
  Eigen::Vector3d contactVel(1.0, 1.0, 1.0);
  
  double contactVelGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  double comVelJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  /*
  double comVxJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  double comVyJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  double comVzJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  */
  
  size_t ii = 0;
  // (1.1) Loop over the caurse grids:
  for(const auto & vel : caurseContactVelocityGrids_)
  {

    contactVel = vel * impactNormal;
    contactVelGrids[ii] = vel;

    updatePiParams_(impactNormal, getParams().eePosition, contactVel);
    // Reset the x component of the impact velocity
    twoDimModelPtr_->updateParams(getPlanarImpactParams());
    // std::cout<<"twoDimModelPtr_->updated params "<<std::endl;
    twoDimModelPtr_->update();

    planarSolutionTo3DPushWall_(velCases_[vel]);

    // Compute the projected com vel jump.
    
    comVelJumpGrids[ii] = impactNormal.transpose() * velCases_[vel].linearVelJump;
    /*
    comVxJumpGrids[ii] = velCases_[vel].linearVelJump.x();
    comVyJumpGrids[ii] = velCases_[vel].linearVelJump.y();
    comVzJumpGrids[ii] = velCases_[vel].linearVelJump.z();
    */

    ii++;
  }

  // (1.2) curve fitting. 
  //gradientApproximationMulti_(contactVelGrids, comVxJumpGrids, comVyJumpGrids, comVzJumpGrids, robotPostImpactStates_.gradient);
  //gradientApproximation_(contactVelGrids, comVxJumpGrids, comVyJumpGrids, comVzJumpGrids, robotPostImpactStates_.gradient);
  // gradientApproximationCalc_(impactNormal, jumpDirection, c1);
  
  gradientApproximation_(contactVelGrids, comVelJumpGrids, robotPostImpactStates_.gradient);

  c1 = robotPostImpactStates_.gradient.c1;

  //jumpDirection = robotPostImpactStates_.gradient.coe;
  double mean = (getTwoDimModelBridgeParams().gradientParams.lowerVelBound + getTwoDimModelBridgeParams().gradientParams.upperVelBound)/2.0;
  jumpDirection = velCases_.equal_range(mean).first->second.linearVelJump.normalized();


}

void TwoDimModelBridge::gradientApproximation_(
		  const double * contactVelGrids, 
		  const double * comVelGrids,
		  fittingParams & params
		  )
{

  double c0, c1, cov11, sumsq;
  // double cov00, cov01;

  // gsl_fit_linear(contactVelGrids, 1, comVelGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

  c0 = 0.0;
  gsl_fit_mul(contactVelGrids, 1, comVelGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &c1, &cov11, &sumsq);
  params.c0 = c0;
  params.c1 = c1;
}


void TwoDimModelBridge::gradientApproximation_(
		  const double * contactVelGrids, 
		  const double * comVxGrids,
		  const double * comVyGrids,
		  const double * comVzGrids,
		  fittingParams & params
		  )
{

  fittingParams cx, cy, cz;

  gsl_fit_mul(contactVelGrids, 1, comVxGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &cx.c1, &cx.cov11, &cx.sumsq);

  gsl_fit_mul(contactVelGrids, 1, comVyGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &cy.c1, &cy.cov11, &cy.sumsq);

  gsl_fit_mul(contactVelGrids, 1, comVzGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &cz.c1, &cz.cov11, &cz.sumsq);

  params.coe.x() = cx.c1;
  params.coe.y() = cy.c1;
  params.coe.z() = cz.c1;

}

void TwoDimModelBridge::gradientApproximationMulti_(
		  const double * contactVelGrids, 
		  const double * comVxGrids,
		  const double * comVyGrids,
		  const double * comVzGrids,
		  fittingParams & params
		  )
{
  size_t i;
  const size_t n(getTwoDimModelBridgeParams().gradientParams.numCaurseGrid);

  const size_t dim = 3; /* linear fit */
  gsl_matrix *ComVelData, *cov;

  // contactVel(i) = vx * coe(0) + vy * coe(1) + vz * coe(2)
  
  //gsl_vector *vx,*vy,*vz;
  gsl_vector *contactVel, *coe;


  ComVelData = gsl_matrix_alloc (n, dim);
  //vx = gsl_vector_alloc (n);
  //vy = gsl_vector_alloc (n);
  //vz = gsl_vector_alloc (n);

  contactVel = gsl_vector_alloc (n);

  coe = gsl_vector_alloc (dim);
  cov = gsl_matrix_alloc (dim, dim);

  
  /* construct design matrix ComVelData for linear fit */
  for (i = 0; i < n; ++i)
  {
    // Set the post-impact COM velocity jump
    gsl_matrix_set (ComVelData, i, 0, comVxGrids[i]);
    gsl_matrix_set (ComVelData, i, 1, comVyGrids[i]);
    gsl_matrix_set (ComVelData, i, 2, comVzGrids[i]);

    // Set the contact velocity
    gsl_vector_set(contactVel, i, contactVelGrids[i]);
  }

  // Allocate the work space
  gsl_multifit_robust_workspace * workSpace
    = gsl_multifit_robust_alloc (gsl_multifit_robust_bisquare, ComVelData->size1, ComVelData->size2);

  // Perform the fitting
  gsl_multifit_robust (ComVelData, contactVel, coe, cov, workSpace);
  gsl_multifit_robust_free (workSpace);

  #define COE(i) (gsl_vector_get(coe,(i)))
  #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  // Extract the coefficients

  //std::cout<<"the fitted params are: "<<COE(0)<<", "<<COE(1)<<", "<<COE(2)<<std::endl;
  params.coe.x() = 1.0/COE(0);
  params.coe.y() = 1.0/COE(1);
  params.coe.z() = 1.0/COE(2);

  //std::cout<<"the gradient is: "<<params.coe.transpose()<<std::endl;

  params.cov(0,0) = COV(0,0);
  params.cov(0,1) = COV(0,1);
  params.cov(0,2) = COV(0,2);

  params.cov(1,0) = COV(1,0);
  params.cov(1,1) = COV(1,1);
  params.cov(1,2) = COV(1,2);

  params.cov(2,0) = COV(2,0);
  params.cov(2,1) = COV(2,1);
  params.cov(2,2) = COV(2,2);

}





Eigen::Vector3d TwoDimModelBridge::projectToGround_(const Eigen::Vector3d & direction)
{
  auto normalizedVec = direction.normalized();
  auto groundNormal = Eigen::Vector3d::UnitZ();
  
  auto result = (groundNormal.cross(normalizedVec)).cross(groundNormal);
  result.normalize();

  return result; 

}
void TwoDimModelBridge::initializeGradientApproximation_()
{

  double caurseStepSize = (getTwoDimModelBridgeParams().gradientParams.upperVelBound
                           - getTwoDimModelBridgeParams().gradientParams.lowerVelBound)
                          / (double)getTwoDimModelBridgeParams().gradientParams.numCaurseGrid;
  std::cout << RoboticsUtils::info <<"The caurse step size is: " << caurseStepSize << RoboticsUtils::reset << std::endl;

  for(double vel = getTwoDimModelBridgeParams().gradientParams.lowerVelBound;
      vel <= getTwoDimModelBridgeParams().gradientParams.upperVelBound; vel += caurseStepSize)
  {
    caurseContactVelocityGrids_.push_back(vel);
    PostImpactStates state;
    velCases_.insert(std::make_pair(vel, state));
    // velCases_.insert(std::make_pair(vel, Eigen::Vector3d::Zero()));

    std::cout << RoboticsUtils::info << "Added contact vel: " << vel << RoboticsUtils::reset << std::endl;
  }
}

void TwoDimModelBridge::printResult()
{
  std::cout << RoboticsUtils::info << "I_nr is: " << twoDimModelPtr_->getSolution().I_nr << RoboticsUtils::reset << std::endl;
  std::cout << RoboticsUtils::info << "I_nc is: " << twoDimModelPtr_->getSolution().I_nc << RoboticsUtils::reset << std::endl;
  std::cout << RoboticsUtils::info << "I_r is: " << twoDimModelPtr_->getSolution().I_r.transpose() <<  RoboticsUtils::reset << std::endl;
  std::cout << RoboticsUtils::info << "The impulse is: " << robotPostImpactStates_.impulse <<  RoboticsUtils::reset << std::endl;

  std::cout << RoboticsUtils::info << "The post-imapct linear velocity is: " << twoDimModelPtr_->getImpactBodies().first.postVel <<  RoboticsUtils::reset << std::endl;
}
void TwoDimModelBridge::printPIParams()
{
  std::cout << RoboticsUtils::info << "The rotation angle is: " << rotationAngle_ << RoboticsUtils::reset << std::endl;
  std::cout << RoboticsUtils::info << "The nu is: " << piParams_.nu.transpose()<< RoboticsUtils::reset  << std::endl;
  std::cout << RoboticsUtils::info << "The tu is: " << piParams_.tu.transpose()<< RoboticsUtils::reset  << std::endl;

  std::cout << RoboticsUtils::info << "The rotated contact point is: " << piParams_.contactPoint.transpose() << RoboticsUtils::reset << std::endl;

  std::cout << RoboticsUtils::info << "The rotated com is: " << piParams_.batParams.com.transpose() << RoboticsUtils::reset << std::endl;

  std::cout << RoboticsUtils::info << "The bat inertia is: " << piParams_.batParams.inertia<< RoboticsUtils::reset  << std::endl;

  std::cout << RoboticsUtils::info << "The bat preimpact vel is: " << piParams_.batParams.preVel << RoboticsUtils::reset << std::endl;

  std::cout << RoboticsUtils::info << "The bat preimpact angular vel is: " << piParams_.batParams.preW << RoboticsUtils::reset << std::endl;
}

void TwoDimModelBridge::updatePiParams_(const Eigen::Vector3d & in,
                                        const Eigen::Vector3d vc,
                                        const Eigen::Vector3d & impactLinearVel)
{
  // (1) Update the normal and tangential unit vectors
  // Compute the angle
  // rotationAngle_= atan2(in.z(), in.y());
  rotationAngle_ = 0.0;
  // std::cout<<green<<"The rotation angle is: "<<rotationAngle_<<std::endl;
  // Update the 2*3 rotation matrix:
  rotation_(0, 0) = 1.0;
  rotation_(1, 1) = cos(rotationAngle_);
  rotationFull_(1, 1) = cos(rotationAngle_);
  rotation_(1, 2) = -sin(rotationAngle_);
  rotationFull_(1, 2) = -sin(rotationAngle_);

  rotationFull_(2, 1) = sin(rotationAngle_);
  rotationFull_(2, 2) = cos(rotationAngle_);

  piParams_.nu = rotation_ * in;
  // Eigen::Vector2d rotatedZ = rotation_*Eigen::Vector3d::UnitZ();
  piParams_.tu(0) = -piParams_.nu(1);
  piParams_.tu(1) = piParams_.nu(0);

  // std::cout<<green<<"The nu is: "<<piParams_.nu.transpose()<<std::endl;
  // std::cout<<green<<"The tu is: "<<piParams_.tu.transpose()<<std::endl;

  // (2) Contact Point:
  piParams_.contactPoint = rotation_ * vc;

  // std::cout<<green<<"The rotated contact point is: "<<piParams_.contactPoint.transpose()<<std::endl;
  // (3) Parmams of the bat and the object:
  switch(getCase_())
  {
    case TwoDimModelCase::PushWall:
      // Update the parameters using the Push-Wall assumptions
      paramUpdatePushWall_(impactLinearVel);
      break;
    default:
      RoboticsUtils::throw_runtime_error("The assumptions are not set for the TwoDimModelBridge.", __FILE__, __LINE__);
  }
}

void TwoDimModelBridge::paramUpdatePushWall_(const Eigen::Vector3d & impactLinearVel)
{
  // Bat is supposed to be the robot:
  // (1) Robot
  piParams_.batParams.com = rotation_ * getRobot()->com();
  // std::cout<<green<<"The rotated com is: "<<piParams_.batParams.com.transpose()<<std::endl;
  piParams_.batParams.mass = getRobot()->mass();
  // Get the z-axis diagonal element:
  piParams_.batParams.inertia = (rotationFull_ * rCentroidalInertia_)(2, 2);
  // std::cout<<green<<"The bat inertia is: "<<piParams_.batParams.inertia<<std::endl;
  
  piParams_.batParams.preVel = rotation_ * impactLinearVel;

  // std::cout<<green<<"The bat preimpact vel is: "<<piParams_.batParams.preVel<<std::endl;
  // Get the z-axis average angular velocity:
  piParams_.batParams.preW = rAverageAngularVel_(2);
  // std::cout<<green<<"The bat preimpact angular vel is: "<<piParams_.batParams.preW<<std::endl;

  piParams_.batParams.name = "robot";

  // Object is suppose to be the wall:

  // (2) Wall
  piParams_.objectParams.preVel << 0.0, 0.0;
  piParams_.objectParams.preW = 0.0;
  piParams_.objectParams.name = "wall";

  // mass and inertia of the wall are set to be infinite.
  // piParams_.objectParams.mass = std::numeric_limits<double>::infinity();
  piParams_.objectParams.mass = std::numeric_limits<double>::max();
  // piParams_.objectParams.inertia = std::numeric_limits<double>::infinity();
  piParams_.objectParams.inertia = std::numeric_limits<double>::max();

  // com of the wall is set to be the contact point such that r = cp - com == 0.
  // Suppose that piParams_.contactPoint is already set.
  piParams_.objectParams.com = piParams_.contactPoint;

  if(getTwoDimModelBridgeParams().debug)
  {
    printPIParams();
  }
}

void TwoDimModelBridge::planarSolutionTo3DPushWall_(PostImpactStates & input)
{
  // Convert the post-impact impulse:
  // The robot applies the impulse "I", thus it receives impulse "-I".
  input.impulse = -rotation_.transpose() * twoDimModelPtr_->getSolution().I_r;
  // std::cout<<"I_nr is: "<< twoDimModelPtr_->getSolution().I_nr<<std::endl;
  // std::cout<<"I_nc is: "<< twoDimModelPtr_->getSolution().I_nc<<std::endl;
  // std::cout<<"I_r is: "<< twoDimModelPtr_->getSolution().I_r.transpose()<<std::endl;
  // std::cout<<"The impulse is: "<< robotPostImpactStates_.impulse<<std::endl;
  // Convert the post-impact velocities:
  // robot:
  input.linearVel = rotation_.transpose() * twoDimModelPtr_->getImpactBodies().first.postVel;
  input.linearVelJump = input.impulse / getRobot()->mass();

  // std::cout<<"The post-imapct linear velocity is: "<< twoDimModelPtr_->getImpactBodies().first.postVel<<std::endl;

  //  Compute: wJump = (1 / getParams().inertia) * cross2(r, I_r);
  Eigen::Vector3d rb;

  rb = getParams().eePosition - getRobot()->com();

  input.anguleVelJump = rCentroidalInertia_.inverse() * rb.cross(input.impulse);

  input.anguleVel = rAverageAngularVel_ + input.anguleVelJump;
}

void TwoDimModelBridge::planarSolutionTo3DPushWall_()
{

  planarSolutionTo3DPushWall_(robotPostImpactStates_);
  /*
  // Convert the post-impact impulse:
  // The robot applies the impulse "I", thus it receives impulse "-I".
  robotPostImpactStates_.impulse = - rotation_.transpose()*twoDimModelPtr_->getSolution().I_r;
  //std::cout<<"I_nr is: "<< twoDimModelPtr_->getSolution().I_nr<<std::endl;
  //std::cout<<"I_nc is: "<< twoDimModelPtr_->getSolution().I_nc<<std::endl;
  //std::cout<<"I_r is: "<< twoDimModelPtr_->getSolution().I_r.transpose()<<std::endl;
  //std::cout<<"The impulse is: "<< robotPostImpactStates_.impulse<<std::endl;
  // Convert the post-impact velocities:
  // robot:
  robotPostImpactStates_.linearVel = rotation_.transpose()*twoDimModelPtr_->getImpactBodies().first.postVel;
  robotPostImpactStates_.linearVelJump = robotPostImpactStates_.impulse/getRobot().mass();

  // std::cout<<"The post-imapct linear velocity is: "<< twoDimModelPtr_->getImpactBodies().first.postVel<<std::endl;

  //  Compute: wJump = (1 / getParams().inertia) * cross2(r, I_r);
  Eigen::Vector3d rb = vcPtr_->getVirtualContactPoint() - vcParams_.com;
  robotPostImpactStates_.anguleVelJump = rCentroidalInertia_.inverse()*rb.cross(robotPostImpactStates_.impulse);

  robotPostImpactStates_.anguleVel = rAverageAngularVel_ + robotPostImpactStates_.anguleVelJump;

  */
  if(getTwoDimModelBridgeParams().debug)
  {
    printResult();
  }
}

void TwoDimModelBridge::positiveCalc_(const double & c1, const double & zmpLowerBoundNorm, const double & zmpUpperBoundNorm, double & maxContactVel, double & minContactVel)
{
  double omega = getRobot()->omega();

  minContactVel = (1.0/c1) * (zmpLowerBoundNorm) * omega;
  maxContactVel = (1.0/c1) * (zmpUpperBoundNorm) * omega;


}

void TwoDimModelBridge::negativeCalc_(const double & c1, const double & zmpLowerBoundNorm, const double & zmpUpperBoundNorm, double & maxContactVel, double & minContactVel)
{
  double omega = getRobot()->omega();

  minContactVel = (1.0/c1) * -(zmpUpperBoundNorm) * omega;
  maxContactVel = (1.0/c1) * (- zmpLowerBoundNorm) * omega;

}

void TwoDimModelBridge::planarSolutionTo3D_()
{
  switch(getCase_())
  {
    case TwoDimModelCase::PushWall:
      // Convert the solution using the Push-Wall assumptions
      planarSolutionTo3DPushWall_();
      break;
    default:
      RoboticsUtils::throw_runtime_error("The assumptions are not set for the TwoDimModelBridge.", __FILE__, __LINE__);
  }
}

const PostImpactStates & TwoDimModelBridge::getObjectPostImpactStates()
{
  switch(getCase_())
  {
    case TwoDimModelCase::PushWall:
      // In this case the object(wall)  is supposed to be stationary.
    RoboticsUtils::throw_runtime_error(
          "In the PushWall case, the wall is stationary. Thus there is no need to check its post-impact states.", __FILE__, __LINE__);
    default:
      return objectPostImpactStates_;
  }
}

const PostImpactStates & ImpactDynamicsModel::getRobotPostImpactStates()
{
  return robotPostImpactStates_;
}

const PostImpactStates & ImpactDynamicsModel::getObjectPostImpactStates()
{
  return objectPostImpactStates_;
}


} // namespace mc_impact 
