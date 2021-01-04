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
#include <Eigen/Dense>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
//# include <Eigen/QR>
#include <mc_control/fsm/Controller.h>
#include <mc_rbdyn/Robots.h>

#include <RBDyn/Momentum.h>

#include <eigen-lssol/LSSOL_QP.h>

#include "ImpactDynamics/ImpactDynamicsModelInterface.h"

#include "ImpactDynamics/OneDimImpactModel.h"
#include "mi_iniEquality.h"
#include "mi_osd.h"
#include "mi_utils.h"
#include "mi_fid_impulse.h"
#include "mi_contactConstraint.h"
#include "mi_unilateralContactConstraint.h"


#include "mi_velIniEquality.h"
#include "mi_velBalance.h"
#include "mi_velInvOsdEquality.h"
#include "mi_velJsdEquality.h"

namespace mc_impact
{

class mi_velEstimator
/*!
 * \brief Estimates the post-impact state jumps.
 */
{
public:
  mi_velEstimator(const mc_rbdyn::Robot & simRobot,
		  const std::shared_ptr<mc_impact::TwoDimModelBridge> twoDimFidModelPtr_,
                 const std::shared_ptr<mi_osd> osdPtr,
		 const std::shared_ptr<mi_qpEstimator> qpEstimator,
                 const struct qpEstimatorParameter params);
  ~mi_velEstimator();
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

  inline const Eigen::MatrixXd & getJacobianTwoDeltaAlpha()
  {
    return jacobianTwoDeltaAlpha_;
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

  /*!
   * \return the Jacobian of the angular momentum derivative jump 
   *
   \\ clang-format off
   \f$ \Delta \dot{L} = \frac{1}{\delta t} \mathcal{J}_{\Delta \dot{L}} \dot{q}_{k+1}  $\f
   \\ clang-format on 
   */
  inline const Eigen::MatrixXd & getJacobianDeltaAM()
  {
    return amJumpJacobian_;
  }
  /*!
   * \return the Jacobian of the liner momentum derivative jump: 
   *
   \\ clang-format off
   \f$ \Delta \dot{P} = \frac{1}{\delta t} \mathcal{J}_{\Delta \dot{P}} \dot{q}_{k+1}  $\f
   \\ clang-format on
   */
  inline const Eigen::MatrixXd & getJacobianDeltaLM()
  {
    return lmJumpJacobian_;
  }

  inline const Eigen::VectorXd & getTauJump() const
  {
    return tauJump_;
  }

  inline const Eigen::VectorXd & getJointVelJump()
  {
    return jointVelJump_;
  }

  /*! \return COM velocity jump \f$ \Delta \dot{c} $\f
   */
  inline const Eigen::Vector3d & getCOMVelJump()
  {
    return comVelJump_;
  }

  /*! \return linear momentum derivative jump \f$ \Delta \dot{P} $\f
   */
  inline const Eigen::Vector3d & getLMJump()
  {
    return lmJump_;
  }

  /*! \return angular momentum derivative jump \f$ \Delta \dot{L} $\f

   */
  inline const Eigen::Vector3d & getAMJump()
  {
    return amJump_;
  }

  const endEffector & getEndeffector(const std::string & name);
  void print() const;
  void print(const std::string & eeName);
  const std::shared_ptr<mi_impactModel> getImpactModel(const std::string & eeName);
  inline const std::map<std::string, std::shared_ptr<mi_impactModel>> & getImpactModels()
  {
    return impactModels_;
  }
  inline const qpEstimatorParameter & getEstimatorParams() const
  {
    return params_;
  }
  inline const std::shared_ptr<mi_osd> getOsd() const
  {
    return osdPtr_;
  }
  inline int getEeNum() const
  {
    return static_cast<int>(endEffectors_.size());
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
    return solverTime_;
  }
  inline double getQweight() const
  {
    return params_.Qweight;
  }
  
  inline void setHostCtl(mc_control::fsm::Controller * ctlPtr)
  {

    if(hostCtlPtr_ == nullptr)
    {
      hostCtlPtr_ = ctlPtr;
    }
    else
    {
      RoboticsUtils::throw_runtime_error("The host fsm controller of the velEstimator: " + getEstimatorParams().name
                               + " is already set!", __FILE__, __LINE__);
    }
  }
  

  /*! \brief Add the GUI entries
   *   Require to set the host fsm controller first
   **/
  void addMcRtcGuiItems();

  /*! \brief Add the log entries
   *   Require to set the host fsm controller first
   **/
  void logImpulseEstimations();

  inline std::shared_ptr<rbd::CentroidalMomentumMatrix> getCmm() const
  {
    return cmmPtr_;
  }
  inline std::shared_ptr<mc_impact::TwoDimModelBridge> getFidModel() const
  {
    return twoDimFidModelPtr_;
  }

private:
  const mc_rbdyn::Robot & simRobot_;

  const std::shared_ptr<mc_impact::TwoDimModelBridge> twoDimFidModelPtr_;
  const std::shared_ptr<mi_osd> osdPtr_;
  const std::shared_ptr<mi_qpEstimator> getQpEstimator_()
  {
    return qpEstimator_;
  }

  const std::shared_ptr<mi_qpEstimator> qpEstimator_;

  endEffector & getEndeffector_(const std::string & name);
  qpEstimatorParameter params_;

  void update_();

  std::vector<std::string> guiEntries_;
  std::vector<std::string> logEntries_;
  void removeImpulseEstimations_();
  void removeMcRtcGuiItems();
  mc_control::fsm::Controller * hostCtlPtr_ = nullptr;

  std::shared_ptr<rbd::CentroidalMomentumMatrix> cmmPtr_;
  inline mc_control::fsm::Controller * getHostCtl_()
  {
    if(hostCtlPtr_ != nullptr)
    {
      return hostCtlPtr_;
    }
    else
    {
      RoboticsUtils::throw_runtime_error("The host fsm controller of the velEstimator: " + getEstimatorParams().name
                               + " is not set!", __FILE__, __LINE__);
    }
  }

  bool osdContactEe_(const std::string & eeName);
  void updateImpactModels_(const std::map<std::string, Eigen::Vector3d> & surfaceNormals);
  void updateImpactModels_();

  std::map<std::string, std::shared_ptr<mi_impactModel>> impactModels_;

  int nameToIndex_(const std::string & eeName);
  bool addEndeffector_(const std::string & eeName, const bool & fromOsd = false);

  void initializeQP_();

  Eigen::VectorXd cl_, cu_;

  Eigen::VectorXd xl_, xu_;

  Eigen::LSSOL_QP solver_;

  Eigen::MatrixXd Q_;
  Eigen::MatrixXd C_;
  Eigen::VectorXd p_;

  std::vector<std::shared_ptr<mi_equality>> eqConstraints_;

  // Add the inequality entries
  std::vector<std::shared_ptr<mi_inEquality>> ieqConstraints_;


  void solveWeightedEqQp_(const Eigen::MatrixXd & Q_,
                          const Eigen::VectorXd & p_,
                          const Eigen::MatrixXd & C_,
                          const Eigen::VectorXd & cu_,
                          Eigen::VectorXd & solution);

  void solveEqQp_(const Eigen::MatrixXd & Q_,
                  const Eigen::VectorXd & p_,
                  const Eigen::MatrixXd & C_,
                  const Eigen::VectorXd & cu_,
                  Eigen::VectorXd & solution);

  void readEeJacobiansSolution_(const Eigen::VectorXd & solutionVariables);

  int numVar_ = 0;
  int numEq_ = 0;
  int numIeq_ = 0;

  inline int getNumVar_() const
  {
    return numVar_;
  }
  inline int getNumEq_() const
  {
    return numEq_;
  }
    inline int getNumIeq_() const
  {
    return numIeq_;
  }
  inline int getNumCon_() const
  {
    return getNumEq_() + getNumIeq_();
  }

  std::map<std::string, endEffector> endEffectors_;
  Eigen::VectorXd jointVelJump_, tauJump_;

  Eigen::Vector3d comVelJump_ = Eigen::Vector3d::Zero();

  Eigen::Vector3d amJump_ = Eigen::Vector3d::Zero();

  ///< amdJump_ = amJumpJacobian_ * jointVelocity_next
  Eigen::MatrixXd amJumpJacobian_;

  Eigen::Vector3d lmJump_ = Eigen::Vector3d::Zero();

  ///< lmdJump_ = lmJumpJacobian_ * jointVelocity_next
  Eigen::MatrixXd lmJumpJacobian_;

  Eigen::MatrixXd jacobianDeltaAlpha_;
  Eigen::MatrixXd jacobianTwoDeltaAlpha_;
  Eigen::MatrixXd jacobianDeltaTau_;
  std::vector<Eigen::MatrixXd> vector_A_dagger_;
  Eigen::MatrixXd tempInv_;

  double solverTime_;
  double structTime_;


};

} // namespace mc_impact
