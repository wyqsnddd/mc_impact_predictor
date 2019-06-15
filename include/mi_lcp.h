# pragma once 

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
  std::vector<double>& solveLCP(const Eigen::MatrixXd & H, const Eigen::VectorXd & d );
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
		int dim
		);

  ~mi_lcp(){}
  void update();
  void update(std::map<std::string, Eigen::Vector3d> contactSurfaceNormals);
  inline const mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }
  protected: 

  void update_(const Eigen::MatrixXd & Jacobian, const Eigen::MatrixXd & JacobianDot);

  Eigen::VectorXd beta_; ///< \beta = \dot{J}\dot{q} - JM^{-1}C
  Eigen::VectorXd d_; ///< d = JM^{-1}\tau + \beta 

  mc_rbdyn::Robot & robot_;
  const std::shared_ptr<mi_osd> & osdPtr_;

  lcp_solver solver_;
  inline const int & getDim_()
  {
   return dim_; 
  }
  int dim_;
   
};
