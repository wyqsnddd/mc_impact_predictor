# pragma once 
#include <iomanip>
# include "mi_osd.h" 
# include <nlopt.hpp>
# include <math.h>
#include <iostream>
#include <Eigen/StdVector>

struct quadraticObjData{
 const Eigen::MatrixXd H;
 const Eigen::VectorXd d;
};

struct lcp_solver{
  static double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *obj_data);
  //static double positiveForceConstraint( unsigned n, const double *x, double *grad, void *data);
  std::vector<double>& solveLCP(const Eigen::MatrixXd & H, const Eigen::VectorXd & d, const std::string & solverName, double convergenceThreshold);
 nlopt::result result;
 std::vector<double> solution;
};

class  mi_lcp
/** \brief  We solve the LCP to predict the contact force for end-effectors with an established contact. 
 */
{
  public: 
  mi_lcp( mc_rbdyn::Robot & robot,
		const std::shared_ptr<mi_osd> & osdPtr,
		int dim,
	  	const std::string & solverName,
		double convergenceThreshold
		);

  ~mi_lcp(){}
  void update();
  void update(std::map<std::string, Eigen::Vector3d> contactSurfaceNormals);
  inline const mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }
  const Eigen::VectorXd  & getPredictedContactForce(const std::string & bodyName);
  inline const int & getDim()
  {
   return dim_; 
  }
  protected: 

  void update_(const Eigen::MatrixXd & Jacobian, const Eigen::MatrixXd & JacobianDot);

  Eigen::VectorXd beta_; ///< \beta = \dot{J}\dot{q} - JM^{-1}C
  Eigen::VectorXd d_; ///< d = JM^{-1}\tau + \beta 

  mc_rbdyn::Robot & robot_;
  const std::shared_ptr<mi_osd> & osdPtr_;

  lcp_solver solver_;
  
  int dim_;
  std::map<std::string, Eigen::VectorXd> predictedContactForce_;
  inline const std::string & getSolver_()
  {
    return solverName_; 
  } 
  std::string solverName_;
  inline const double & getThreshold_()
  {
    return convergenceThreshold_; 
  }
  double convergenceThreshold_;
};
