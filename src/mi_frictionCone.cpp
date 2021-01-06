#include "mi_frictionCone.h"

namespace mc_impact
{

mi_frictionCone::mi_frictionCone(const std::shared_ptr<mi_osd> osdPtr) : mi_inEquality(osdPtr)
{

  reset_();

  std::cout << RoboticsUtils::alarm << "Initialized the fulfilling friction cone constraint" << RoboticsUtils::reset
            << std::endl;
}

void mi_frictionCone::reset_()
{
  // Use the endEffectors with contacts
  int dof = getOsd_()->getDof();
  int nContacts = static_cast<int>(getOsd_()->getContactEes().size());
  int nEe = static_cast<int>(getOsd_()->getEeNum());
  // int nContactEe = static_cast<int>(getOsd_()->getContactNum());
  // int nImpactEe = static_cast<int>(impactModels_.size());

  // We need to adapt to the possible dim of the Jacobian:
  // 1: directional Jacobian
  // 3: linear Jacobian
  // 6: full Jacobian
  int dim = 3;
  int rowNum = 5;

  // The dimension of A: is dim \times established contacts
  A_.resize(rowNum * nContacts, dof + dim * (nEe));
  A_.setZero();

  // b_ stays zero
  b_.resize(rowNum * nContacts);
  b_.setZero();

  // Assume z-axis is the contact normal direction.
  Eigen::MatrixXd matrixA;
  matrixA.resize(rowNum, dim);
  matrixA.setZero();

  // Ix - 1/miu*Iz <= 0
  matrixA(0, 0) = 1;
  matrixA(0, 2) = -1 / getMiu_();
  // -Ix - 1/miu*Iz <= 0
  matrixA(1, 0) = -1;
  matrixA(1, 2) = -1 / getMiu_();
  // Iy - 1/miu*Iz <= 0
  matrixA(2, 1) = 1;
  matrixA(2, 2) = -1 / getMiu_();
  matrixA(3, 1) = -1;
  matrixA(3, 2) = -1 / getMiu_();
  // Iz >= 0
  // matrixA(4, 2) = 1;

  int count(0);

  for(auto & contact : getOsd_()->getContactEes())
  {
    int eeIndex = getOsd_()->nameToIndex_(contact);
    A_.block(count * rowNum, dof + eeIndex * dim, rowNum, dim) = matrixA;
    count++;
  }
}

void mi_frictionCone::update() {}

} // namespace mc_impact
