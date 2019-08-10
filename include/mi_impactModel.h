# pragma once

# include <iostream>
# include <Eigen/StdVector>

# include "mi_osd.h"

class mi_impactModel
{
  public: 
   mi_impactModel(
		   const mc_rbdyn::Robot & simRobot,
		   const std::shared_ptr<mi_osd> & osdPtr,
		   std::string iBodyName="r_wrist",
		   double iDuration = 0.005,
		   double timeStep = 0.005,
		   double coeF = 0.2,
		   double coeR = 0.8,
		   int dim = 3
		 ): simRobot_(simRobot), osdPtr_(osdPtr), impactBodyName_(iBodyName), impactDuration_(iDuration), timeStep_(timeStep), coeFrictionDeduction_(coeF), coeRes_(coeR), dim_(dim)
   {
   }

   ~mi_impactModel(){}

  inline int  getDim() const
  {
   return dim_; 
  }
  inline const std::string & getImpactBody() const
  {
   return impactBodyName_; 
  }
  inline double getTimeStep() const
  {
    return timeStep_;
  }
  inline double getImpactDuration() const
  {
    return impactDuration_;
  }
  inline double getCoeRes() const
  {
    return coeRes_;
  }
  inline double getCoeFricDe() const
  {
    return coeFrictionDeduction_;
  }

  inline const Eigen::VectorXd & getEeVelocityJump() const
  {
    return deltaV_;
  }
  inline const Eigen::MatrixXd & getProjector() const 
  {
    return reductionProjector_;
  }
  void update(const Eigen::Vector3d & surfaceNormal); 

  private:
  const mc_rbdyn::Robot & simRobot_;
  const std::shared_ptr<mi_osd> & osdPtr_;
  std::string impactBodyName_;
  double impactDuration_;
  double timeStep_;
  double coeFrictionDeduction_;
  double coeRes_;
  int dim_;

  Eigen::VectorXd deltaV_ = Eigen::VectorXd::Zero(3);
  Eigen::MatrixXd reductionProjector_ = Eigen::MatrixXd::Zero(3,3);

};
