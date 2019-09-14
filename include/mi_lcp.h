#pragma once
#include "mi_osd.h"
#include <Eigen/StdVector>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <nlopt.hpp>

namespace mc_impact
{

struct quadraticObjData
{
  const Eigen::MatrixXd H;
  const Eigen::VectorXd d;
};

struct lcp_solver
{
  static double objFunction(const std::vector<double> & x, std::vector<double> & grad, void * obj_data);
  // static double positiveForceConstraint( unsigned n, const double *x, double *grad, void *data);
  std::vector<double> & solveLCP(const Eigen::MatrixXd & H,
                                 const Eigen::VectorXd & d,
                                 const std::string & solverName,
                                 double convergenceThreshold);
  nlopt::result result;
  std::vector<double> solution;
};

class mi_lcp
/** \brief  We solve the LCP to predict the contact force for end-effectors with an established contact.
 */
{
public:
  mi_lcp(const mc_rbdyn::Robot & simRobot,
         const mc_rbdyn::Robot & realRobot,
         const std::shared_ptr<mi_osd> & osdPtr,
         int dim,
         const std::string & solverName,
         double convergenceThreshold);

  ~mi_lcp() {}
  void update();
  void update(const std::map<std::string, Eigen::Vector3d> & contactSurfaceNormals);
  inline const mc_rbdyn::Robot & getRealRobot()
  {
    return realRobot_;
  }
  inline const mc_rbdyn::Robot & getSimRobot()
  {
    return simRobot_;
  }
  const Eigen::VectorXd & getPredictedContactForce(const std::string & bodyName);
  inline int getDim()
  {
    return dim_;
  }

protected:
  void update_(const Eigen::MatrixXd & Jacobian, const Eigen::MatrixXd & JacobianDot);

  Eigen::VectorXd beta_; ///< \beta = \dot{J}\dot{q} - JM^{-1}C
  Eigen::VectorXd d_; ///< d = JM^{-1}\tau + \beta

  const mc_rbdyn::Robot & simRobot_;
  const mc_rbdyn::Robot & realRobot_;
  const std::shared_ptr<mi_osd> osdPtr_;

  lcp_solver solver_;

  int dim_;
  std::map<std::string, Eigen::VectorXd> predictedContactForce_;
  inline const std::string & getSolver_()
  {
    return solverName_;
  }
  std::string solverName_;
  inline double getThreshold_()
  {
    return convergenceThreshold_;
  }
  double convergenceThreshold_;
};
} // namespace mc_impact
