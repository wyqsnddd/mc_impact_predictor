#pragma once

#include <mc_control/fsm/Controller.h>

#include <RBDyn/Momentum.h>

#include <TwoDimModel/TwoDimModel.h>
#include <VirtualContactPoint/SolveSemiAxes.h>
#include <VirtualContactPoint/VirtualContactPoint.h>

// Header file for GNU Scientific Library: Least squares fit.
#include <gsl/gsl_fit.h>

#include "mi_utils.h"

namespace mc_impact
{

struct PostImpactStates
{
  Eigen::Vector3d linearVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d linearVelJump = Eigen::Vector3d::Zero();
  Eigen::Vector3d anguleVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d anguleVelJump = Eigen::Vector3d::Zero();
  Eigen::Vector3d impulse = Eigen::Vector3d::Zero();
  double c0 = 0.0;
  double c1 = 0.0;
};

class ImpactDynamicsModel
{
public:
  ImpactDynamicsModel(const mc_rbdyn::Robot & simRobot, const ImpactModelParams & params);

  ~ImpactDynamicsModel() {}

  virtual const PostImpactStates & getRobotPostImpactStates();
  virtual const PostImpactStates & getObjectPostImpactStates();

  inline const mc_rbdyn::Robot & getRobot()
  {
    return simRobot_;
  }

  inline const ImpactModelParams & getParams()
  {
    return params_;
  }

  virtual void update(const Eigen::Vector3d & impactNormal, const Eigen::Vector3d & impactLinearVel) = 0;
  // virtual void update() = 0;
protected:
  const mc_rbdyn::Robot & simRobot_;
  ImpactModelParams params_;
  PostImpactStates robotPostImpactStates_;
  PostImpactStates objectPostImpactStates_;
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
  bool useVirtualContact = true;
  bool useComVel = true;
  bool debug = false;
  bool gradientApproximation = true;
  GradientApproximationParams gradientParams;
  ImpactModelParams modelParams;
};

class TwoDimModelBridge : public ImpactDynamicsModel
/*! \brief We always assume the inertial frame
 * In this case, we assume k
 */
{
public:
  TwoDimModelBridge(const mc_rbdyn::Robot & simRobot, const TwoDimModelBridgeParams & params);

  ~TwoDimModelBridge() {}
  const PostImpactStates & getObjectPostImpactStates() override;

  void update(const Eigen::Vector3d & impactNormal, const Eigen::Vector3d & impactLinearVel) override;
  // void update() override;

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
    return params_;
  }
  inline void setHostCtl(mc_control::fsm::Controller * ctlPtr)
  {

    if(hostCtlPtr_ == nullptr)
    {
      hostCtlPtr_ = ctlPtr;
    }
    else
    {
      throw std::runtime_error("The host fsm controller is already set!");
    }
  }
  void logImpulseEstimations();

  void removeImpulseEstimations();
  const Eigen::Vector3d & getAverageLinerVel()
  {
    return rAverageLinearVel_;
  }
  const Eigen::Vector3d & getAverageAngularVel()
  {
    return rAverageAngularVel_;
  }
  void printVcParams();
  void printPIParams();
  void printResult();


  /*! \return the moment of inertia of the entire robot in its centroidal frame.
   */
  const Eigen::Matrix3d & getRobotCentroidalInertia()
  {
    return rCentroidalInertia_; 
  }

protected:
  TwoDimModelBridgeParams params_;
  // bool useVirtualContact_;

  double rotationAngle_;

  std::vector<std::string> logEntries_;
  mc_control::fsm::Controller * hostCtlPtr_ = nullptr;
  inline mc_control::fsm::Controller * getHostCtl_()
  {
    if(hostCtlPtr_ != nullptr)
    {
      return hostCtlPtr_;
    }
    else
    {
      throw std::runtime_error("The host fsm controller is not set!");
    }
  }
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

  std::shared_ptr<FIDynamics::VirtualContactPoint> vcPtr_;
  FIDynamics::VcParams vcParams_;
  inline const FIDynamics::VcParams & getVcParams_()
  {
    return vcParams_;
  }

  std::shared_ptr<FIDynamics::SolveSemiAxes> ssaPtr_;

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
  std::vector<double> caurseContactVelocityGrids_;
  // std::vector<double> fineContactVelocityGrids_;
  void initializeGradientApproximation_();

}; // end of the TwoDimModelBridge

/*
class OneDimModelBridge : public ImpactDynamicsModel
{
public:
OneDimModelBridge(const mc_rbdyn::Robot & simRobot,
    const ImpactModelParams & params);

void update() override;

protected:

};
*/

} // namespace mc_impact
