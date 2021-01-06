/* Copyright 2019 CNRS-UM LIRMM
 *
 * \author Yuquan Wang, Arnaud Tanguy
 *
 *
 *
 * mc_impact_predictor is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * mc_impact_predictor is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with mc_impact_predictor. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "mi_osd.h"
#include "mi_utils.h"
#include <Eigen/StdVector>
#include <RoboticsUtils/utils.h>
#include <chrono>
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

struct LcpSolver
{
  static double objFunction(const std::vector<double> & x, std::vector<double> & grad, void * obj_data);
  // static double positiveForceConstraint( unsigned n, const double *x, double *grad, void *data);
  std::vector<double> & solveLCP(const Eigen::MatrixXd & H,
                                 const Eigen::VectorXd & d,
                                 const std::string & solverName,
                                 double convergenceThreshold);
  nlopt::result result;
  std::vector<double> solution;
  double solverTime; ///< Time to solve the LP in milliseconds.
};

class mi_lcp
/** \brief  We solve the LCP to predict the contact force for end-effectors with an established contact.
 */
{
public:
  mi_lcp(const mc_rbdyn::Robot & simRobot,
         const mc_rbdyn::Robot & realRobot,
         const std::shared_ptr<mi_osd> osdPtr,
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
  /*! \brief Time to construct the building blocks in each iteration.
   * \return time in microseconds.
   */
  inline double structTime() const
  {
    return structTime_;
  }

  /*! \brief Time to solve the optimization problem.
   * \return time in microseconds.
   */
  inline double solverTime() const
  {
    return solver_.solverTime;
  }

protected:
  void update_(const Eigen::MatrixXd & Jacobian, const Eigen::MatrixXd & JacobianDot);

  Eigen::VectorXd beta_; ///< \beta = \dot{J}\dot{q} - JM^{-1}C
  Eigen::VectorXd d_; ///< d = JM^{-1}\tau + \beta

  const mc_rbdyn::Robot & simRobot_;
  const mc_rbdyn::Robot & realRobot_;
  const std::shared_ptr<mi_osd> osdPtr_;

  LcpSolver solver_;

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

  double structTime_;
};
} // namespace mc_impact
