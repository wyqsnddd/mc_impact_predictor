# include "mi_balance.h"

namespace mc_impact{

mi_balance::mi_balance(
		const std::shared_ptr<mi_osd> osdPtr,
		const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels
		): mi_equality(osdPtr), impactModels_(impactModels)
{

  reset_();

  std::cout<<red<<"Initialized the impulse balance constraint"<<reset<<std::endl;
}


void mi_balance::reset_()
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

  // The dimension of A is 6 \times dof CMM matrix, in conjunction with nEe blocks 
  A_.resize(6, dof + dim* (nEe));
  A_.setZero();

  // b_ will stay zero
  b_.resize(6);
  b_.setZero();


  // Initialize the centroidal momentum matrix: 
  cmmPtr_ = std::make_shared<rbd::CentroidalMomentumMatrix>(getOsd_()->getRobot().mb());

}

Eigen::Matrix3d mi_balance::crossMatrix(const Eigen::Vector3d & input)
{

  Eigen::Matrix3d skewSymmetricMatrix = Eigen::Matrix3d::Zero();

  skewSymmetricMatrix(0, 1) = -input(2);
  skewSymmetricMatrix(1, 0) = input(2);

  skewSymmetricMatrix(0, 2) = input(1);
  skewSymmetricMatrix(2, 0) = -input(1);

  skewSymmetricMatrix(1, 2) = -input(0);
  skewSymmetricMatrix(2, 1) = input(0);

  return skewSymmetricMatrix;
}

Eigen::MatrixXd mi_balance::forceGraspMatrix(const std::string eeName, const Eigen::Vector3d & reference)
{
   Eigen::MatrixXd graspMatrx;
   graspMatrx.resize(6, 3);
   graspMatrx.setZero();

  auto X_0_c = getOsd_()->getRobot().bodyPosW(eeName);

  auto rotationTranspose = X_0_c.rotation();

  auto translation = X_0_c.translation() - reference;

  graspMatrx.block<3, 3>(0, 0) = rotationTranspose; 

  graspMatrx.block<3, 3>(3, 0) = -rotationTranspose * crossMatrix(translation);

  return graspMatrx;

}

void mi_balance::update()
{

  int dof = getOsd_()->getDof();
  int dim = getOsd_()->getJacobianDim();

  // (0) Fill the centroidal momentum matrix
  
  cmmPtr_ ->computeMatrix(getOsd_()->getRobot().mb(), getOsd_()->getRobot().mbc(), getOsd_()->getRobot().com());

  const Eigen::MatrixXd &  cmmMatrix = cmmPtr_->matrix();


  A_.block(0, 0, 6, dof) = cmmMatrix; 

    // (1) Fill the contact:

  //for(auto idx = getOsd_()->getContactEes().begin(); idx != getOsd_()->getContactEes().end(); ++idx)
  for(auto & contact: getOsd_()->getContactEes())
  {
    int eeIndex = getOsd_()->nameToIndex_(contact);
    // int eeIndex = nameToIndex_(*idx);
    const Eigen::MatrixXd & eeGraspMatrix = forceGraspMatrix(contact, getOsd_()->getRobot().com());
    A_.block(0, dof + eeIndex * dim, 6, dim) = eeGraspMatrix; 
  }

  // (2) Fill the impact bodies: 

  //for(auto idx = impactModels_.begin(); idx != impactModels_.end(); ++idx)
  for(auto & impactModel : impactModels_) 
  {
    int eeIndex = getOsd_()->nameToIndex_(impactModel.first);

    const Eigen::MatrixXd & eeGraspMatrix = forceGraspMatrix(impactModel.first, getOsd_()->getRobot().com());
    A_.block(0, dof + eeIndex* dim, 6, dim) =  eeGraspMatrix;
  }
}







}

