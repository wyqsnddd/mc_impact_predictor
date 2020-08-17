#include <TwoDimModel/TwoDimModel.h>
#include <VirtualContactPoint/VirtualContactPoint.h>
#include <VirtualContactPoint/SolveSemiAxes.h>

#include <RBDyn/Momentum.h>
#include "mi_utils.h"

namespace mc_impact 
{


class ImpactDynamicsModel
{
public:

ImpactDynamicsModel(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params);

virtual const Eigen::VectorXd & getPostImpactVel();

inline const mc_rbdyn::Robot & getRobot()
{
  return simRobot_;
}

inline const ImpactModelParams  & getParams()
{
  return params_;
}

virtual void update(const Eigen::Vector3d & impactNormal); 
virtual void update();
protected:

const mc_rbdyn::Robot & simRobot_;
ImpactModelParams params_;

};


class TwoDimModelBridge : public ImpactDynamicsModel
/*! \brief We always assume the inertial frame
 */
{
public:
TwoDimModelBridge(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params);

const Eigen::VectorXd & getPostImpactVel() override;

void update(const Eigen::Vector3d & impactNormal) override;
void update() override;

const FIDynamics::PIParams & getPlanarImpactParams()
{
  return piParams_;
}

protected:

void updatePlanarImpactParams_();
FIDynamics::PIParams piParams_;

std::shared_ptr<FIDynamics::TwoDimModel> twoDimModelPtr_;

std::shared_ptr<FIDynamics::VirtualContactPoint> vcPtr_;
FIDynamics::VcParams vcParams_;
inline const FIDynamics::VcParams & getVcParams_()
{
  return vcParams_;
}

std::shared_ptr<FIDynamics::SolveSemiAxes> ssaPtr_;

};


class OneDimModelBridge : public ImpactDynamicsModel
{
public:
OneDimModelBridge(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params);

const Eigen::VectorXd & getPostImpactVel() override;

void update() override;

protected:

};



} // End of NameSpace
