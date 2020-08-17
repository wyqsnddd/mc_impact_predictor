#include <TwoDimModel/TwoDimModel.h>
#include <VirtualContactPoint/VirtualContactPoint.h>
#include <VirtualContactPoint/SolveSemiAxes.h>
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
protected:

const mc_rbdyn::Robot & simRobot_;
ImpactModelParams params_;

};


class TwoDimModelBridge : public ImpactDynamicsModel
{
public:
TwoDimModelBridge(const mc_rbdyn::Robot & simRobot,
		const ImpactModelParams & params);

const Eigen::VectorXd & getPostImpactVel() override;

void update();
protected:

FIDynamics::PIParams params_;

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

protected:

};



} // End of NameSpace
