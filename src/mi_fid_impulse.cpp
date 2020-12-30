#include "mi_fid_impulse.h"

namespace mc_impact
{

mi_fid_impulse::mi_fid_impulse(const std::shared_ptr<mi_osd> osdPtr,
		const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels,

		const std::shared_ptr<rbd::CentroidalMomentumMatrix> cmmPtr,
		const std::shared_ptr<mc_impact::TwoDimModelBridge> twoDimModel) : mi_equality(osdPtr), impactModels_(impactModels), cmmPtr_(cmmPtr), twoDimFidModelPtr_(twoDimModel) 
{

  reset_();

  std::cout << red << "Initialized the fid impulse constraint" << reset << std::endl;
}

void mi_fid_impulse::reset_()
{
  // Use the endEffectors with contacts
  int dof = getOsd_()->getDof();
  int nEe = static_cast<int>(getOsd_()->getEeNum());
  // int nContactEe = static_cast<int>(getOsd_()->getContactNum());
  // int nImpactEe = static_cast<int>(impactModels_.size());

  // We need to adapt to the possible dim of the Jacobian:
  // 1: directional Jacobian
  // 3: linear Jacobian
  // 6: full Jacobian
  int dim = getOsd_()->getJacobianDim();

  // The dimension of A: is 6 \times dof + dim * (Ee blocks)
  A_.resize(6, dof + dim * (nEe));
  A_.setZero();

  // b_ collects the impulse 
  b_.resize(6);
  b_.setZero();

  // Initialize the centroidal momentum matrix:
  cmmPtr_ = std::make_shared<rbd::CentroidalMomentumMatrix>(getOsd_()->getRobot().mb());

}

void mi_fid_impulse::update()
{

  int dof = getOsd_()->getDof();
  int dim = getOsd_()->getJacobianDim();

  // (0) Fill the impulse 
  // Linear Momentum Jump: mass * linear vel jump
  b_.segment<3> (3) = getOsd_()->getRobot().mass() * getFidModel()->getRobotPostImpactStates().linearVelJump;

  // Angular Momentum Jump: inertia * angular vel jump
  b_.segment<3> (0) = getFidModel()->getRobotCentroidalInertia() * getFidModel()->getRobotPostImpactStates().anguleVelJump;


  cmmPtr_->sComputeMatrix(getOsd_()->getRobot().mb(), getOsd_()->getRobot().mbc(), getOsd_()->getRobot().com());


  //const Eigen::MatrixXd & cmmMatrix = cmmPtr_->matrix();

  //A_.block(0, 0, 6, dof) = cmmMatrix;

  // (1) Fill the contact:

  
  
  // for(auto idx = getOsd_()->getContactEes().begin(); idx != getOsd_()->getContactEes().end(); ++idx)
  for(auto & contact : getOsd_()->getContactEes())
  {
    int eeIndex = getOsd_()->nameToIndex_(contact);
    // int eeIndex = nameToIndex_(*idx);
    const Eigen::MatrixXd & eeGraspMatrix = getOsd_()->forceGraspMatrix(contact, getOsd_()->getRobot().com());
    A_.block(0, dof + eeIndex * dim, 6, dim) = eeGraspMatrix;
  }

  // (2) Fill the impact bodies:

  // for(auto idx = impactModels_.begin(); idx != impactModels_.end(); ++idx)
  for(auto & impactModel : impactModels_)
  {
    int eeIndex = getOsd_()->nameToIndex_(impactModel.first);

    const Eigen::MatrixXd & eeGraspMatrix = getOsd_()->forceGraspMatrix(impactModel.first, getOsd_()->getRobot().com());
    A_.block(0, dof + eeIndex * dim, 6, dim) = eeGraspMatrix;
  }
  
  
}

} // namespace mc_impact
