#pragma once

#include <TwoDimModel/TwoDimModel.h>
#include <McDynamicStability/McZMPArea.h>
#include <RobotInterface/RobotInterface.h>

// Header file for GNU Scientific Library: Least squares fit.
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>

#include "../mi_utils.h"
#include "utils.h"

namespace mc_impact 
{

enum class TwoDimModelCase
{
  PushWall
};
struct fittingParams
{
  Eigen::Vector3d coe = Eigen::Vector3d::Zero();
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();

  double c0 = 0.0;
  double c1 = 0.0;
  
  double cov00 = 0.0;
  double cov01 = 0.0;
  double cov11 = 0.0;

  double sumsq; 
};

struct PostImpactStates
{
  Eigen::Vector3d linearVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d linearVelJump = Eigen::Vector3d::Zero();
  Eigen::Vector3d anguleVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d anguleVelJump = Eigen::Vector3d::Zero();
  Eigen::Vector3d impulse = Eigen::Vector3d::Zero();
  fittingParams gradient;
  double c = 0.0;
  //double c0 = 0.0;
  //double c1 = 0.0;
};

struct GradientApproximationParams
{
  double upperVelBound = 0.8; // Meter/Second
  double lowerVelBound = 0.0; // Meter/Second
  size_t numCaurseGrid = 20; // Number of grids between upper and lower Velocity Bound.

  // Compute for the points v +- (m * size)
  // m = velBound/stepSize
  double velBound = 0.1; // Meter/Second
  double stepSize = 0.01; // Meter/Second
};

struct TwoDimModelBridgeParams
{
  std::string name = "TwoDimModel";
  //bool useVirtualContact = true;
  //bool useComVel = true;
  bool debug = false;
  bool gradientApproximation = true;
  GradientApproximationParams gradientParams;
  //ImpactModelParams modelParams;
};

struct StandingStabilityBounds
{
  double maxContactVel = 0.0;
  double minContactVel = 0.0;

  double maxCOMVel = 0.0;
  double minCOMVel = 0.0;

};


struct StandingStabilityParams
{
  Eigen::Vector2d  zmpLowerBound = Eigen::Vector2d::Zero();
  Eigen::Vector2d  zmpUpperBound = Eigen::Vector2d::Zero();

  //double zmpLowerBoundNorm = -0.1; 
  //double zmpUpperBoundNorm= 0.1; 

  ///< Direction of the post-impact com-Velocity jump
  Eigen::Vector3d jumpDirection = Eigen::Vector3d::Zero();

  double pseudoCOMVel = 0.05; 

  double c1 = 0.0;
  double omega = 0.0;

  double pOne = -1.0;
  double pTwo = -0.2;

  StandingStabilityBounds weakBounds;
  StandingStabilityBounds strongBounds;

};


class ImpactDynamicsModel
{
public:
  ImpactDynamicsModel(const std::shared_ptr<RobotInterface::Robot> robotPtr,
 const ImpactModelParams & params);

  virtual ~ImpactDynamicsModel() {
    std::cout<<RobotInterface::info<<"Destructing ImpactDynamicsModel." <<RobotInterface::reset<<std::endl;
  }

  virtual const PostImpactStates & getRobotPostImpactStates();
  virtual const PostImpactStates & getObjectPostImpactStates();

  inline const std::shared_ptr<RobotInterface::Robot> getRobot()
  {
    return robotPtr_;
  }

  inline const ImpactModelParams & getParams()
  {
    return params_;
  }

  virtual void update(const Eigen::Vector3d & impactNormal, const Eigen::Vector3d & impactLinearVel) = 0;
  virtual void update(const Eigen::Vector3d & impactNormal, const std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_, StandingStabilityParams & params) = 0;
  // virtual void update() = 0;
  
  /*! \brief Set the friction coefficient. 
   */
  virtual void setFrictionCoefficient(const double & miu) = 0;

  /*! \brief Set the restitution coefficient. 
   */
  virtual void setRestitutionCoefficient(const double & e) = 0;

  /*! \brief Set the name of the impact body / surface
   */
  virtual void setImpactSurfaceName(const std::string & iBody) = 0;


protected:
  std::shared_ptr<RobotInterface::Robot> robotPtr_;
  ImpactModelParams params_;
  PostImpactStates robotPostImpactStates_;
  PostImpactStates objectPostImpactStates_;
};

class TwoDimModelBridge : public ImpactDynamicsModel
/*! \brief We always assume the inertial frame
 * In this case, we assume k
 */
{
public:
  TwoDimModelBridge(
		  const std::shared_ptr<RobotInterface::Robot> robotPtr,
                  const ImpactModelParams & params,
		  const TwoDimModelBridgeParams & brigeParams);

  virtual ~TwoDimModelBridge() {
    std::cout<<RobotInterface::info<<"Destructing TwoDimModelBridge." <<RobotInterface::reset<<std::endl;
  }
  const PostImpactStates & getObjectPostImpactStates() override;

  void update(const Eigen::Vector3d & impactNormal, const Eigen::Vector3d & impactLinearVel) override;

  void update(const Eigen::Vector3d & impactNormal, const std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_, StandingStabilityParams & params) override;

  /*
  inline const Eigen::Vector3d & getImpulse()
  {
   return postImpactImpulse_;
  }
  */
  const FIDynamics::PIParams & getPlanarImpactParams()
  {
    return piParams_;
  }
  inline const TwoDimModelBridgeParams & getTwoDimModelBridgeParams()
  {
    return bridgeParams_;
  }
   
  

  const Eigen::Vector3d & getAverageLinerVel()
  {
    return rAverageLinearVel_;
  }
  const Eigen::Vector3d & getAverageAngularVel()
  {
    return rAverageAngularVel_;
  }
  void printPIParams();
  void printResult();


  /*! \return the moment of inertia of the entire robot in its centroidal frame.
   */
  const Eigen::Matrix3d & getRobotCentroidalInertia()
  {
    return rCentroidalInertia_; 
  }

  void setFrictionCoefficient(const double & miu) override
  {
    //twoDimModelPtr_->setFrictionCoefficient(miu);
    piParams_.miu = miu;
  }

  /*! \brief Set the energetic restitution coefficient. 
   */
  void setRestitutionCoefficient(const double & e) override
  {
    //twoDimModelPtr_->setEnergeticRestitutionCoefficient(e);
    piParams_.e = e;
  }

  void setImpactSurfaceName(const std::string & iSurface) override
  {
    params_.iSurfaceName= iSurface;
  }
protected:
  TwoDimModelBridgeParams bridgeParams_;

  double rotationAngle_;

  std::vector<std::string> logEntries_;

  // Compute the planar impact parameters using 3D data.
  /*!
   * \param in The impact normal
   * \param vc The virtual contact point
   */
  void updatePiParams_(const Eigen::Vector3d & in, const Eigen::Vector3d vc, const Eigen::Vector3d & impactLinearVel);
  FIDynamics::PIParams piParams_;

  Eigen::Matrix<double, 2, 3> rotation_;
  Eigen::Matrix<double, 3, 3> rotationFull_;
  // Eigen::Vector3d postImpactImpulse_ = Eigen::Vector3d::Zero();

  Eigen::Matrix3d rCentroidalInertia_; ///< Robot centroidal inertia.
  Eigen::Vector3d rAverageAngularVel_; ///< Robot average angular velocity
  Eigen::Vector3d rAverageLinearVel_; ///< Robot average angular velocity


  std::shared_ptr<FIDynamics::TwoDimModel> twoDimModelPtr_;

  // Convert the 2D solution to 3D
  void planarSolutionTo3D_();

  /*! \brief indicates what assumptions are used.
   */
  inline const TwoDimModelCase & getCase_()
  {
    return case_;
  }

  TwoDimModelCase case_ = TwoDimModelCase::PushWall;
  /*! \brief update the piParmas with the Push-Wall assumption: object mass and moment of inertia are infinite.
   */
  void paramUpdatePushWall_(const Eigen::Vector3d & impactLinearVel);
  void planarSolutionTo3DPushWall_();

  void planarSolutionTo3DPushWall_(PostImpactStates & input);

  // First: impact velocity, Second: Corresponding PostImpactStates
  // std::map<double, Eigen::Vector3d> velCases_;
  std::map<double, PostImpactStates> velCases_;

  void gradientApproximationMulti_(
		  const double * contactVelGrids, 
		  const double * comVxGrids,
		  const double * comVyGrids,
		  const double * comVzGrids,
		  fittingParams & params
		  );


  void gradientApproximation_(
		  const double * contactVelGrids, 
		  const double * comVxGrids,
		  const double * comVyGrids,
		  const double * comVzGrids,
		  fittingParams & params
		  );

  std::vector<double> caurseContactVelocityGrids_;

  // std::vector<double> fineContactVelocityGrids_;
  void initializeGradientApproximation_();

  /*! \brief Project the 3D impact normal on the ground plane perpendicular to [0,0,1]
   */
  Eigen::Vector3d projectToGround_(const Eigen::Vector3d & direction);

  void negativeCalc_(const double & c1, const double & zmpLowerBoundNorm, const double & zmpUpperBoundNorm, double & maxContactVel, double & minContactVel);
  void positiveCalc_(const double & c1, const double & zmpLowerBoundNorm, const double & zmpUpperBoundNorm, double & maxContactVel, double & minContactVel);

  void velSaturation_(double & inputVel);
  /*! \brief Compute the gradient of the post-impact COM velocity jump w.r.t. the contact velocity
 */
  void computeGradient_(const Eigen::Vector3d & impactNormal, StandingStabilityParams & params);

/*! \brief Computes the maximum feasible contact vel along the given impact normal.  The ZMP bound along the impact normal is also known.
 */
  void computeMaxContactVel_(const Eigen::Vector3d & impactNormal, const Eigen::Vector2d & zmpLowerBound, const Eigen::Vector2d & zmpUpperBound, StandingStabilityParams & params);


}; // end of the TwoDimModelBridge

} // namespace mc_impact
