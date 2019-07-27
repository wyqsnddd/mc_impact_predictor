# pragma once 
#include <iomanip>
# include "mi_osd.h" 
# include <nlopt.hpp>
# include <math.h>
#include <iostream>
#include <Eigen/StdVector>



#include "mi_osd.h"
#include <mc_rbdyn/Robots.h>

//struct constraintObjData{
// const Eigen::MatrixXd H;
// const Eigen::VectorXd d;
//};
class mi_qpEstimator;

struct qpEstimatorParameter{
  std::string impactBodyName = "r_wrist";
  double impactDuration;
  double coeFrictionDeduction = 0.2;
  double coeRes = 0.8;
  int dim = 3;
  std::string solverName = "nlopt::LD_CCSAQ";
  double convergenceThreshold = 0.01;
};


struct quadraticObjData{
 const Eigen::MatrixXd H;
};

struct endEffector{
 //std::string name;
 //unsigned dim; // dimension 
 size_t startingNumber; // Starting number in the vector
 Eigen::VectorXd estimatedImpulse;
};

struct qp_solver{
  static void jsdImpulseConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data);
  static void osdImpulseConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data);
  static void iniConstraintFunction( unsigned constraintDim, double *result, unsigned stateDim, const double * x, double* grad, void * f_data);
  static double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *obj_data);
  std::vector<double>& solveQP( const mi_qpEstimator* estimatorPtr);
 nlopt::result result;
 std::vector<double> solution;
};

class mi_qpEstimator{
  public: 
  mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> & osdPtr,
		const struct qpEstimatorParameter params
		//const std::string & impactBodyName,
		//int  dim,
	  	//const std::string & solverName,
		//double convergenceThreshold
		);
  ~mi_qpEstimator(){}
  void update(const Eigen::Vector3d & surfaceNormal);
  inline const mc_rbdyn::Robot & getSimRobot()
  {
    return simRobot_;
  }
  
  const Eigen::VectorXd  & getPredictedImpulse(const std::string & bodyName) const;

  
  void addEndeffector(std::string eeName);

  inline const int & getDim() const
  {
   return params_.dim; 
  }
  inline const int & getDof() const
  {
   return getOsd_()->getDof(); 
  }
  inline const std::string & getImpactBody() const
  {
   return params_.impactBodyName; 
  }
  inline const Eigen::VectorXd & getEeVelocityJump() const
  {
    return deltaV_;
  }
  private: 
  const mc_rbdyn::Robot & simRobot_;
  const std::shared_ptr<mi_osd> & osdPtr_;
  inline const std::shared_ptr<mi_osd> & getOsd_() const
  {
    return osdPtr_;
  }
  qpEstimatorParameter params_;

  Eigen::VectorXd deltaV_; ///< End-effector velocity jump \f$ \delta \v_k \f$
  /*
  std::string impactBodyName_;
  double coeFrictionDeduction_;
  double coeRes_;
  int dim_;
*/

  qp_solver solver_;
  std::map<std::string, endEffector> endEffectors_;

  //void addOptVariables_( const std::string& name, const unsigned & dim);

  const endEffector & getEndeffector_( const std::string& name) const;
  inline const double & getImpactDuration_() const
  {
    return params_.impactDuration;
  }
  inline const double & getCoeRes_() const
  {
    return params_.coeRes;
  }
  inline const double & getCoeFricDe_() const
  {
    return params_.coeFrictionDeduction;
  }

  inline const std::string & getSolver_() const
  {
    return params_.solverName; 
  }
  inline const double & getThreshold_() const
  {
    return params_.convergenceThreshold; 
  }
  friend struct qp_solver; 
  friend struct impactDynamicsConstraintData; 
};
