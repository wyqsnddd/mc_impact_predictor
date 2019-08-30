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
#include <Eigen/StdVector>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
//# include <Eigen/QR>
#include <mc_rbdyn/Robots.h>

#include "mi_impactModel.h"
#include "mi_iniEquality.h"
#include "mi_invOsdEquality.h"
#include "mi_jsdEquality.h"
#include "mi_osd.h"
#include "mi_utils.h"
#include <Eigen/Dense>
#include <eigen-lssol/LSSOL_QP.h>

namespace mc_impact
{

class mi_qpEstimator
{
public:
  mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
                 const std::shared_ptr<mi_osd> osdPtr,
                 const struct qpEstimatorParameter params);
  ~mi_qpEstimator() {}
  void update(const std::map<std::string, Eigen::Vector3d> & surfaceNormals);
  void update();

  inline const mc_rbdyn::Robot & getSimRobot()
  {
    return simRobot_;
  }

  inline int getDof() const
  {
    return getOsd()->getDof();
  }

  inline double getQweight() const
  {
    return params_.Qweight;
  }

  inline const Eigen::MatrixXd & getJacobianDeltaAlpha()
  {
    return jacobianDeltaAlpha_;
  }
  inline const Eigen::MatrixXd & getJacobianDeltaTau()
  {
    return jacobianDeltaTau_;
  }

  inline const Eigen::MatrixXd & getJacobianDeltaF(const std::string & eeName)
  {
    return getEndeffector(eeName).jacobianDeltaF;
  }
  inline const Eigen::VectorXd & getTauJump() const
  {
    return tauJump_;
  }
  inline const Eigen::VectorXd & getJointVelJump()
  {
    return jointVelJump_;
  }

  const endEffector & getEndeffector(const std::string & name);
  void print() const;
  void print(const std::string & eeName);
  const std::shared_ptr<mi_impactModel> & getImpactModel(const std::string & eeName);
  inline const std::map<std::string, std::shared_ptr<mi_impactModel>> & getImpactModels()
  {
    return impactModels_;
  }
  inline const qpEstimatorParameter & getEstimatorParams()
  {
    return params_;
  }
  inline const std::shared_ptr<mi_osd> & getOsd() const
  {
    return osdPtr_;
  }
  inline int getEeNum()
  {
    return static_cast<int>(endEffectors_.size());
  }

private:
  const mc_rbdyn::Robot & simRobot_;
  const std::shared_ptr<mi_osd> osdPtr_;
  endEffector & getEndeffector_(const std::string & name);
  qpEstimatorParameter params_;
  void update_();

  bool osdContactEe_(const std::string & eeName);
  void updateImpactModels_(const std::map<std::string, Eigen::Vector3d> & surfaceNormals);
  void updateImpactModels_();

  std::map<std::string, std::shared_ptr<mi_impactModel>> impactModels_;

  int nameToIndex_(const std::string & eeName);
  bool addEndeffector_(const std::string & eeName, const bool & fromOsd = false);

  void initializeQP_();

  Eigen::MatrixXd Q_;

  Eigen::MatrixXd C_;

  Eigen::VectorXd cl_, cu_;

  Eigen::VectorXd p_;

  Eigen::VectorXd xl_, xu_;
  Eigen::LSSOL_QP solver_;

  std::vector<std::shared_ptr<mi_equality>> eqConstraints_;

  void solveEqQp_(const Eigen::MatrixXd & Q_,
                  const Eigen::VectorXd & p_,
                  const Eigen::MatrixXd & C_,
                  const Eigen::VectorXd & cu_,
                  Eigen::VectorXd & solution);

  inline int getNumVar_() const
  {
    return numVar_;
  }
  int numVar_;

  inline int getNumEq_() const
  {
    return numEq_;
  }
  int numEq_;

  std::map<std::string, endEffector> endEffectors_;
  Eigen::VectorXd jointVelJump_, tauJump_;

  Eigen::MatrixXd jacobianDeltaAlpha_;
  Eigen::MatrixXd jacobianDeltaTau_;
  std::vector<Eigen::MatrixXd> vector_A_dagger_;
  Eigen::MatrixXd tempInv_;
};
} // namespace mc_impact
