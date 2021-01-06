#include "mi_velBalance.h"

namespace mc_impact
{

mi_velBalance::mi_velBalance(const std::shared_ptr<mi_osd> osdPtr,
                             const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels,
                             const std::shared_ptr<rbd::CentroidalMomentumMatrix> cmmPtr,
                             const std::shared_ptr<mc_impact::TwoDimModelBridge> twoDimmodel)
: mi_equality(osdPtr), impactModels_(impactModels), cmmPtr_(cmmPtr), twoDimFidModelPtr_(twoDimmodel)
{

  reset_();

  std::cout << RoboticsUtils::alarm << "Initialized the momentum balance constraint" << RoboticsUtils::reset
            << std::endl;
}

void mi_velBalance::reset_()
{
  // Use the endEffectors with contacts
  int dof = getOsd_()->getDof();

  // The dimension of A is 6 \times dof CMM matrix, in conjunction with nEe blocks
  A_.resize(6, dof);
  A_.setZero();

  // b_ will stay zero
  b_.resize(6);
  b_.setZero();

  // Initialize the centroidal momentum matrix:
  cmmPtr_ = std::make_shared<rbd::CentroidalMomentumMatrix>(getOsd_()->getRobot().mb());
}

void mi_velBalance::update()
{

  // (0) Fill the impulse
  // Linear Momentum Jump: mass * linear vel jump
  b_.segment<3>(3) = getOsd_()->getRobot().mass() * getFidModel()->getRobotPostImpactStates().linearVelJump;

  // Angular Momentum Jump: inertia * angular vel jump
  b_.segment<3>(0) =
      getFidModel()->getRobotCentroidalInertia() * getFidModel()->getRobotPostImpactStates().anguleVelJump;

  // (0) Fill the centroidal momentum matrix

  cmmPtr_->sComputeMatrix(getOsd_()->getRobot().mb(), getOsd_()->getRobot().mbc(), getOsd_()->getRobot().com());

  const Eigen::MatrixXd & cmmMatrix = cmmPtr_->matrix();

  A_ = cmmMatrix;
  // std::cout<<red<<"The CMM matrix is: "<<std::endl<< cyan<<cmmMatrix<<reset<<std::endl;
}

} // namespace mc_impact
