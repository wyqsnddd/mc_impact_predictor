#pragma once

#include <RobotInterface/RobotInterface.h>
#include <RoboticsUtils/utils.h>
#include <TwoDimModel/TwoDimModel.h>

// Header file for GNU Scientific Library: Least squares fit.
#include "../mi_utils.h"
#include "utils.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>

namespace mc_impact
{

enum class TwoDimModelCase
{
  PushWall
};
struct fittingParams
{
  // Eigen::Vector3d coe = Eigen::Vector3d::Zero();
  //Eigen::Matrix3d cov = Eigen::Matrix3d::Identity();

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
    //double c = 0.0;
  // double c0 = 0.0;
  // double c1 = 0.0;
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

  fittingParams gradientX;
  fittingParams gradientY;

};

struct TwoDimModelBridgeParams
{
  std::string name = "TwoDimModel";
  // bool useVirtualContact = true;
  // bool useComVel = true;
  bool debug = false;
  bool gradientApproximation = true;
  GradientApproximationParams gradientParams;
  // ImpactModelParams modelParams;
};

class ImpactDynamicsModel
{
public:
  ImpactDynamicsModel(const std::shared_ptr<RobotInterface::Robot> robotPtr, const ImpactModelParams & params);

  virtual ~ImpactDynamicsModel()
  {
    std::cout << RoboticsUtils::info << "Destructing ImpactDynamicsModel." << RoboticsUtils::reset << std::endl;
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
  // virtual void update(const Eigen::Vector3d & impactNormal, const
  // std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>> mcZMPAreaPtr_, StandingStabilityParams & params) = 0;
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
  TwoDimModelBridge(const std::shared_ptr<RobotInterface::Robot> robotPtr,
                    const ImpactModelParams & params,
                    const TwoDimModelBridgeParams & brigeParams);

  virtual ~TwoDimModelBridge()
  {
    std::cout << RoboticsUtils::info << "Destructing TwoDimModelBridge." << RoboticsUtils::reset << std::endl;
  }
  const PostImpactStates & getObjectPostImpactStates() override;

  void update(const Eigen::Vector3d & impactNormal, const Eigen::Vector3d & impactLinearVel) override;

  // void update(const Eigen::Vector3d & impactNormal, const std::shared_ptr<mc_impact::McZMPArea<Eigen::Vector2d>>
  // mcZMPAreaPtr_, StandingStabilityParams & params) override;

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
    // twoDimModelPtr_->setFrictionCoefficient(miu);
    piParams_.miu = miu;
  }

  /*! \brief Set the energetic restitution coefficient.
   */
  void setRestitutionCoefficient(const double & e) override
  {
    // twoDimModelPtr_->setEnergeticRestitutionCoefficient(e);
    piParams_.e = e;
  }

  void setImpactSurfaceName(const std::string & iSurface) override
  {
    params_.iSurfaceName = iSurface;
  }

  /*! \brief Compute the gradient of the post-impact COM velocity jump w.r.t. the contact velocity: c * comVelocity =
   * contactVelocity
   */
  void computeGradient(const Eigen::Vector3d & impactNormal, Eigen::Vector3d & jumpDirection);

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

  /*
  void gradientApproximationMulti_(const double * contactVelGrids,
                                   const double * comVxGrids,
                                   const double * comVyGrids,
				   const double * comVzGrids,
                                   fittingParams & params);
				   */

  void gradientApproximation_(const double * contactVelGrids,
                              const double * comVxGrids,
                              const double * comVyGrids,
                              GradientApproximationParams & params);

  //void gradientApproximationCalc_(const Eigen::Vector3d & impactNormal, double & c1);

  //void gradientApproximation_(const double * contactVelGrids, const double * comVGrids, fittingParams & params);

  std::vector<double> caurseContactVelocityGrids_;

  // std::vector<double> fineContactVelocityGrids_;
  void initializeGradientApproximation_();

  /*! \brief Project the 3D impact normal on the ground plane perpendicular to [0,0,1]
   */
  Eigen::Vector3d projectToGround_(const Eigen::Vector3d & direction);

  void negativeCalc_(const double & c1,
                     const double & zmpLowerBoundNorm,
                     const double & zmpUpperBoundNorm,
                     double & maxContactVel,
                     double & minContactVel);
  void positiveCalc_(const double & c1,
                     const double & zmpLowerBoundNorm,
                     const double & zmpUpperBoundNorm,
                     double & maxContactVel,
                     double & minContactVel);

}; // end of the TwoDimModelBridge

} // namespace mc_impact
