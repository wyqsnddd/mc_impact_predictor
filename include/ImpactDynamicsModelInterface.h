#include <TwoDimModel/TwoDimModel.h>
#include <VirtualContactPoint/VirtualContactPoint.h>
#include <VirtualContactPoint/SolveSemiAxes.h>

#include <RBDyn/Momentum.h>
#include "mi_utils.h"

namespace mc_impact 
{

struct PostImpactStates
{
  Eigen::Vector3d linearVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d anguleVel = Eigen::Vector3d::Zero();
  Eigen::Vector3d impulse = Eigen::Vector3d::Zero();
};

class ImpactDynamicsModel
{
public:

ImpactDynamicsModel(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params);

~ImpactDynamicsModel(){}

virtual const PostImpactStates & getRobotPostImpactStates();
virtual const PostImpactStates & getObjectPostImpactStates();

inline const mc_rbdyn::Robot & getRobot()
{
  return simRobot_;
}

inline const ImpactModelParams  & getParams()
{
  return params_;
}

virtual void update(const Eigen::Vector3d & impactNormal) = 0;
virtual void update() = 0;
protected:

const mc_rbdyn::Robot & simRobot_;
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
TwoDimModelBridge(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params);

~TwoDimModelBridge(){}
const PostImpactStates & getObjectPostImpactStates() override;

void update(const Eigen::Vector3d & impactNormal) override;
void update() override;

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

protected:
// Compute the planar impact parameters using 3D data.
/*!
 * \param in The impact normal
 * \param vc The virtual contact point
 */
void updatePiParams_(const Eigen::Vector3d & in, const Eigen::Vector3d vc);
FIDynamics::PIParams piParams_;

Eigen::Matrix<double, 2, 3> rotation_;
Eigen::Matrix<double, 3, 3> rotationFull_;
//Eigen::Vector3d postImpactImpulse_ = Eigen::Vector3d::Zero();

Eigen::Matrix3d rCentroidalInertia_; ///< Robot centroidal inertia.
Eigen::Vector3d rAverageAngularVel_; ///< Robot average angular velocity 

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
void paramUpdatePushWall_();
void planarSolutionTo3DPushWall_();
};


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


} // End of NameSpace
