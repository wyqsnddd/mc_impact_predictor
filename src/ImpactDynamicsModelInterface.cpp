#include "ImpactDynamics/ImpactDynamicsModelInterface.h" 

namespace mc_impact 
{


int sgn(double val) {
    return (double(0) < val) - (val < double(0));
}

ImpactDynamicsModel::ImpactDynamicsModel(
		const std::shared_ptr<RobotInterface::Robot> robotPtr,
		const ImpactModelParams & params)
: robotPtr_(robotPtr), params_(params)
{
  std::cout << RobotInterface::info << "ImpactDynamicsModel is created." << RobotInterface::reset << std::endl;
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

  std::cout << RobotInterface::info << "TwoDimModelBridge is created." << RobotInterface::reset << std::endl;
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
    throw_runtime_error("", __FILE__, __LINE__) ;
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
    /*
    gsl_fit_linear(contactVelGrids, 1, comVelJumpGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		   */

    c0 = 0.0;
    gsl_fit_mul(contactVelGrids, 1, comVelJumpGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &c1, &cov11, &sumsq);

    robotPostImpactStates_.c = c1;
  }
}

void TwoDimModelBridge::update(const Eigen::Vector3d & impactNormal, const std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_, StandingStabilityParams & params)
{
  computeGradient_(impactNormal, params);

  auto cvj_normal = getRobotPostImpactStates().linearVelJump.normalized();

  mcZMPAreaPtr_->computeZMPBound(-cvj_normal, params.zmpLowerBound, params.zmpUpperBound);

  params.zmpUpperDistance = params.zmpUpperBound.norm();
  // If the upper and lower bound are in the same quadrant, then they have the same positive sign, otherwise the lower bound has the negative sign.
   int sameQuadrant = [](const Eigen::Vector2d & lower, const Eigen::Vector2d & upper)->int{  
  bool sameX = sgn(upper.x()) == sgn(lower.x());

  bool sameY = sgn(upper.y()) == sgn(lower.y());

  if (sameX && sameY)
  {
    return 1;
  }
  else 
  {
    return -1;
  }
  }(params.zmpLowerBound, params.zmpUpperBound);

  params.zmpLowerDistance = sameQuadrant * params.zmpLowerBound.norm();

  computeMaxContactVel_(impactNormal, params.zmpLowerDistance, params.zmpUpperDistance, params);
}

void TwoDimModelBridge::computeGradient_(const Eigen::Vector3d & impactNormal, StandingStabilityParams & params)
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
  sva::PTransformd X_0_ee = getRobot()->bodyPosW(getParams().iSurfaceName);

  params_.eePosition = X_0_ee.translation();

  // (1) Approximate the gradient. 
  
  Eigen::Vector3d contactVel(1.0, 1.0, 1.0);
  

  
  double contactVelGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  double comVxJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  double comVyJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  double comVzJumpGrids[getTwoDimModelBridgeParams().gradientParams.numCaurseGrid];
  
  /*
  std::vector<double>contactVelGrids, comVxJumpGrids, comVyJumpGrids, comVzJumpGrids;
  contactVelGrids.resize(getTwoDimModelBridgeParams().gradientParams.numCaurseGrid);
  comVxJumpGrids.resize(getTwoDimModelBridgeParams().gradientParams.numCaurseGrid);
  comVyJumpGrids.resize(getTwoDimModelBridgeParams().gradientParams.numCaurseGrid);
  comVzJumpGrids.resize(getTwoDimModelBridgeParams().gradientParams.numCaurseGrid);
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
    
    //comVelJumpGrids[ii] = impactNormal.transpose() * velCases_[vel].linearVelJump;
    comVxJumpGrids[ii] = velCases_[vel].linearVelJump.x();
    comVyJumpGrids[ii] = velCases_[vel].linearVelJump.y();
    comVzJumpGrids[ii] = velCases_[vel].linearVelJump.z();

    ii++;
  }

  // (1.2) curve fitting. 
  //gradientApproximationMulti_(contactVelGrids, comVxJumpGrids, comVyJumpGrids, comVzJumpGrids, robotPostImpactStates_.gradient);
  gradientApproximation_(contactVelGrids, comVxJumpGrids, comVyJumpGrids, comVzJumpGrids, robotPostImpactStates_.gradient);

  /*
  fittingParams cx, cy, cz;

  gsl_fit_mul(contactVelGrids, 1, comVxJumpGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &cx.c1, &cx.cov11, &cx.sumsq);

  gsl_fit_mul(contactVelGrids, 1, comVyJumpGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &cy.c1, &cy.cov11, &cy.sumsq);

  gsl_fit_mul(contactVelGrids, 1, comVzJumpGrids, 1, getTwoDimModelBridgeParams().gradientParams.numCaurseGrid,
                   &cz.c1, &cz.cov11, &cz.sumsq);
		   */

  // c1 should be negative.
  
  params.jumpDirection = robotPostImpactStates_.gradient.coe;
 // params.jumpDirection.x() = cx.c1;
 // params.jumpDirection.y() = cy.c1;
 // params.jumpDirection.z() = cz.c1;

  params.c1 = params.jumpDirection.norm(); 
  // The sign depends on the angle between the "contact vel" and the "post-impact com velocity jump"
  
 double signTest =  impactNormal.transpose()*params.jumpDirection;

  //std::cout<<"Interface: norm: "<<params.c1<<std::endl;
  params.c1 = -sgn(signTest) * params.c1;
  //std::cout<<"c1 is:: "<<params.c1<<std::endl;
  robotPostImpactStates_.c = params.c1;
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



void TwoDimModelBridge::computeMaxContactVel_(const Eigen::Vector3d & impactNormal, const double & zmpLowerDistance, const double & zmpUpperDistance, StandingStabilityParams & params)
{
  // Compute the strong stability bounds: 
  params.omega = getRobot()->omega();

  params.strongBounds.maxCOMVel = (zmpUpperDistance) * params.omega;
  params.strongBounds.minCOMVel = (zmpLowerDistance) * params.omega;

  auto v1 = (1.0/params.c1) * (params.strongBounds.maxCOMVel - params.pseudoCOMVel);
  auto v2 = (1.0/params.c1) * (params.strongBounds.minCOMVel - params.pseudoCOMVel); 

  params.strongBounds.maxContactVel = v1>v2?v1:v2; 
  velSaturation_(params.strongBounds.maxContactVel);
  
  params.strongBounds.minContactVel = v1>v2?v2:v1; 

  if(params.strongBounds.minContactVel < 0.0)
  {
    params.strongBounds.minContactVel = 0.0;
  }else
  {
    velSaturation_(params.strongBounds.minContactVel);
  }
  

  // Compute the weak stability bounds: 
  /*
  params.weakBounds.maxCOMVel = -(params.zmpUpperBoundNorm/(params.pTwo  + params.pOne)) * params.omega;
  params.weakBounds.minCOMVel = -(params.zmpLowerBoundNorm/(params.pTwo + params.pOne)) * params.omega;

  auto w_v1 = (1.0/c1) * (params.weakBounds.maxCOMVel - params.pseudoCOMVel);
  auto w_v2 = (1.0/c1) * (params.weakBounds.minCOMVel - params.pseudoCOMVel); 

  params.weakBounds.maxContactVel = w_v1>w_v2?w_v1:w_v2; 
  params.weakBounds.minContactVel = w_v1>w_v2?w_v2:w_v1; 
  */

  auto pGain = 1 + params.pOne * params.pTwo;
  auto dGain = -(params.pOne+ params.pTwo) / params.omega;

  //auto projectedCOM = static_cast<double>(impactNormal.transpose() * projectToGround_(getRobot()->com()));
  params.weakBounds.maxCOMVel = (zmpUpperDistance) / dGain;
  params.weakBounds.minCOMVel = (zmpLowerDistance) / dGain;

  auto w_v1 = (1.0/params.c1) * (params.weakBounds.maxCOMVel - params.pseudoCOMVel);
  auto w_v2 = (1.0/params.c1) * (params.weakBounds.minCOMVel - params.pseudoCOMVel); 

  params.weakBounds.maxContactVel = w_v1>w_v2?w_v1:w_v2; 
  velSaturation_(params.weakBounds.maxContactVel);

  params.weakBounds.minContactVel = w_v1>w_v2?w_v2:w_v1; 
  if(params.weakBounds.minContactVel < 0.0)
  {
    params.weakBounds.minContactVel = 0.0;
  }else
  {
    velSaturation_(params.weakBounds.minContactVel);
  }
  /*
  if(c1<0)
  {
    negativeCalc_(c1, zmpLowerBoundNorm, zmpUpperBoundNorm, maxContactVel, minContactVel);
  }
  else if(c1>0)
  {
    positiveCalc_(c1, zmpLowerBoundNorm, zmpUpperBoundNorm, maxContactVel, minContactVel);
  }
  else
  {
    throw std::runtime_error("c1 is zero!!!");
  }
  */
 
  
  if(params.strongBounds.maxContactVel<=0.0 || params.weakBounds.maxContactVel<=0.0)
  {
    throw_runtime_error("maxContactVel is negative!!!", __FILE__, __LINE__);
  }
  //std::cout<<"twodim: c0 is: "<<c0<< std::endl;
  //std::cout<<"twodim: c1 is: "<<c1<< std::endl;
  // (2) Compute the maximum contact vel using the ZMP bound and the approximated gradient.
  // Note that we need  to perform the computation in 2D.
  //double omega = getRobot()->omega();

  //Eigen::Vector3d zmpLowerBound3d( 0.1, 0.0, 0.0);
  //Eigen::Vector3d zmpLowerBound3d( zmpUpperBound(0), zmpUpperBound(1), 0.0);

  //double zmpUpperBoundNorm = zmpUpperBound.norm();
  //double zmpLowerBoundNorm = zmpLowerBound.norm();

  //auto projectedImpactNormal = projectToGround_(impactNormal);
  //std::cout<<"twodim: zmpBound is: "<< zmpLowerBound3d.transpose() << std::endl;

  //std::cout<<"twodim: impact normal is: "<< impactNormal.transpose() << std::endl;
  //std::cout<<"twodim: projected impact normal is: "<< projectedImpactNormal.transpose() << std::endl;
  //std::cout<<"impactNormal.transpose() * ( zmpLowerBound3d + (getRobot()->com() * omega)): "<< projectedImpactNormal.transpose() * ( zmpLowerBound3d + (getRobot()->com() * omega)) <<std::endl;
  //std::cout<<"impactNormal.transpose() * zmpLowerBound3d: "<< projectedImpactNormal.transpose() * zmpLowerBound3d <<std::endl;

  //std::cout<<"impactNormal.transpose() * (getRobot()->com() * omega): "<< projectedImpactNormal.transpose() * ( getRobot()->com() * omega) <<std::endl;

  //maxContactVel = (1.0/c1) * impactNormal.transpose() * ( zmpLowerBound3d + (getRobot()->com() * omega));
  //maxContactVel = (1.0/c1) * projectedImpactNormal.transpose() * (zmpLowerBound3d );
  //maxContactVel = (1.0/c1) * (zmpLowerBound3d.norm()) * omega;
  
  //maxContactVel = (1.0/c1) * (zmpUpperBoundNorm) * omega;
  //minContactVel = (1.0/c1) * -(zmpLowerBoundNorm) * omega;

  
  /*
  if(maxContactVel < 0.0)
  {
    //Eigen::Vector3d zmpUpperBound3d( zmpUpperBound(0), zmpUpperBound(1), 0.0);
    Eigen::Vector3d zmpUpperBound3d( zmpLowerBound(0), zmpLowerBound(1), 0.0);

    //maxContactVel = (1.0/c1) * impactNormal.transpose() * ( zmpUpperBound3d + (getRobot()->com() * omega));
    //maxContactVel = (1.0/c1) * projectedImpactNormal.transpose() * (zmpUpperBound3d );
    
    maxContactVel = (1.0/c1) * (zmpLowerBound3d.norm()) * omega;
  }
  */
  

  //std::cout<<"twodim: maxContactvel is: "<< maxContactVel<< std::endl;
  // (3) Call the two dim model again to fill in the post-impact states 
  
  //std::cout<<"twodim: check point five"<< std::endl;

  
  // We need the fabs to impose maxContactVel the same direction as the impact normal
  
  updatePiParams_(impactNormal, getParams().eePosition, impactNormal * params.weakBounds.maxContactVel);

  twoDimModelPtr_->updateParams(getPlanarImpactParams());
  twoDimModelPtr_->update();
  planarSolutionTo3D_();
  


}

void TwoDimModelBridge::velSaturation_(double & inputVel)
{
 /*
 if(inputVel < 0.0)
 {
   throw_runtime_error("TwoDimModelBridge::The input vel is negative!", __FILE__, __LINE__);
 }
 */

 if(inputVel >= getTwoDimModelBridgeParams().gradientParams.upperVelBound)
 {
   inputVel = getTwoDimModelBridgeParams().gradientParams.upperVelBound;
 }else if(inputVel <= getTwoDimModelBridgeParams().gradientParams.lowerVelBound)
 {
   inputVel = getTwoDimModelBridgeParams().gradientParams.lowerVelBound;
 }


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
  std::cout << RobotInterface::info <<"The caurse step size is: " << caurseStepSize << RobotInterface::reset << std::endl;

  for(double vel = getTwoDimModelBridgeParams().gradientParams.lowerVelBound;
      vel <= getTwoDimModelBridgeParams().gradientParams.upperVelBound; vel += caurseStepSize)
  {
    caurseContactVelocityGrids_.push_back(vel);
    PostImpactStates state;
    velCases_.insert(std::make_pair(vel, state));
    // velCases_.insert(std::make_pair(vel, Eigen::Vector3d::Zero()));

    std::cout << RobotInterface::info << "Added contact vel: " << vel << reset << std::endl;
  }
}

void TwoDimModelBridge::printResult()
{
  std::cout << RobotInterface::info << "I_nr is: " << twoDimModelPtr_->getSolution().I_nr << RobotInterface::reset << std::endl;
  std::cout << RobotInterface::info << "I_nc is: " << twoDimModelPtr_->getSolution().I_nc << RobotInterface::reset << std::endl;
  std::cout << RobotInterface::info << "I_r is: " << twoDimModelPtr_->getSolution().I_r.transpose() <<  RobotInterface::reset << std::endl;
  std::cout << RobotInterface::info << "The impulse is: " << robotPostImpactStates_.impulse <<  RobotInterface::reset << std::endl;

  std::cout << RobotInterface::info << "The post-imapct linear velocity is: " << twoDimModelPtr_->getImpactBodies().first.postVel <<  RobotInterface::reset << std::endl;
}
void TwoDimModelBridge::printPIParams()
{
  std::cout << RobotInterface::info << "The rotation angle is: " << rotationAngle_ << RobotInterface::reset << std::endl;
  std::cout << RobotInterface::info << "The nu is: " << piParams_.nu.transpose()<< RobotInterface::reset  << std::endl;
  std::cout << RobotInterface::info << "The tu is: " << piParams_.tu.transpose()<< RobotInterface::reset  << std::endl;

  std::cout << RobotInterface::info << "The rotated contact point is: " << piParams_.contactPoint.transpose() << RobotInterface::reset << std::endl;

  std::cout << RobotInterface::info << "The rotated com is: " << piParams_.batParams.com.transpose() << RobotInterface::reset << std::endl;

  std::cout << RobotInterface::info << "The bat inertia is: " << piParams_.batParams.inertia<< RobotInterface::reset  << std::endl;

  std::cout << RobotInterface::info << "The bat preimpact vel is: " << piParams_.batParams.preVel << RobotInterface::reset << std::endl;

  std::cout << RobotInterface::info << "The bat preimpact angular vel is: " << piParams_.batParams.preW << RobotInterface::reset << std::endl;
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
      throw_runtime_error("The assumptions are not set for the TwoDimModelBridge.", __FILE__, __LINE__);
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
      throw_runtime_error("The assumptions are not set for the TwoDimModelBridge.", __FILE__, __LINE__);
  }
}

const PostImpactStates & TwoDimModelBridge::getObjectPostImpactStates()
{
  switch(getCase_())
  {
    case TwoDimModelCase::PushWall:
      // In this case the object(wall)  is supposed to be stationary.
      throw_runtime_error(
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


} // namespace RobotInterface 
