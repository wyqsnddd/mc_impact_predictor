# include "mi_jsdEquality.h" 

mi_jsdEquality::mi_jsdEquality(
		const std::shared_ptr<mi_osd> & osdPtr,
		const std::vector<std::string> & impactBodies
		) :  mi_equality(osdPtr), impactBodyNames_(impactBodies)
{
  reset_();
}

void mi_jsdEquality::reset_()
{
  int dof = getOsd_()->getDof();
  int nEe = static_cast<int>(getOsd_()->getEeNum());
  int dim = getOsd_()->getJacobianDim(); 

  A_.resize(dof, dof + nEe*dim);
  A_.setZero();

  // b_ will stay zero
  b_.resize(dof);
  b_.setZero();
}

void mi_jsdEquality::update()
{
  
  int nRow = getOsd_()->getDof();
  int dim = getOsd_()->getJacobianDim(); 
  // (0) Fill the mass matrix
  A_.block(0, 0, nRow, nRow) = getOsd_()->getMassMatrix();

  // (1) Fill the contact bodies
  
  for (auto idx = getOsd_()->getContactEes().begin(); idx !=  getOsd_()->getContactEes().end(); ++idx)
  {
    int eeIndex = getOsd_()->nameToIndex_(*idx);
    A_.block(0, nRow + eeIndex*dim, nRow, dim) = - getOsd_()->getJacobian(*idx).transpose();
  }
  
  // (2) Fill the impact bodies
  for (auto idx = impactBodyNames_.begin(); idx!= impactBodyNames_.end(); ++idx  )
  {
     int eeIndex = getOsd_()->nameToIndex_(*idx);
     A_.block(0, nRow + eeIndex*dim, nRow, dim) = - getOsd_()->getJacobian(*idx).transpose();
  }
 }
