# pragma once 
# include <iomanip>
# include "mi_osd.h" 
# include <math.h>
# include <iostream>
# include <Eigen/StdVector>

# include <mc_rbdyn/Robots.h>
# include <eigen-lssol/LSSOL_QP.h>
# include "mi_jsdEquality.h"
# include "mi_osdEquality.h"
# include "mi_invOsdEquality.h"
# include "mi_iniEquality.h"
# include "mi_impactModel.h"
# include <limits>
//# include <Eigen/QR>
# include <Eigen/Dense>

struct qpEstimatorParameter{
  double Qweight = 20;
  //std::vector<std::string> impactBodyNames={"r_wrist"};
  std::map<std::string, Eigen::Vector3d> impactNameAndNormals;
  double impactDuration = 0.005;
  double timeStep = 0.005;
  double coeFrictionDeduction = 0.2;
  double coeRes = 0.8;
  int dim = 3;
  bool useLagrangeMultiplier = false;
  bool useJsd = true;
  bool useOsd = true;
};

struct endEffector{
 //size_t startingNumber; // Starting number in the vector
 Eigen::Vector3d eeVJump; 
 Eigen::Vector3d estimatedImpulse;
 Eigen::Vector3d estimatedAverageImpulsiveForce;
 Eigen::Vector3d checkForce;
 Eigen::MatrixXd jacobianDeltaF;
};

class mi_qpEstimator{
  public: 
  mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> osdPtr,
		const struct qpEstimatorParameter params
		);
  ~mi_qpEstimator(){
 
  }
  void update(const std::map<std::string, Eigen::Vector3d> & surfaceNormals);
  void update();

  inline const mc_rbdyn::Robot & getSimRobot()
  {
    return simRobot_;
  }
  
  void addEndeffector(std::string eeName);

 
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

  const endEffector & getEndeffector( const std::string& name);
  void print() const;
  void print(const std::string & eeName);
  const std::shared_ptr<mi_impactModel> & getImpactModel(const std::string & eeName);
  inline const std::map<std::string, std::shared_ptr<mi_impactModel> > & getImpactModels()
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
  endEffector & getEndeffector_( const std::string& name);
  qpEstimatorParameter params_;
  void update_();

  void updateImpactModels_(const std::map<std::string, Eigen::Vector3d> & surfaceNormals);
  void updateImpactModels_();
  std::map<std::string, std::shared_ptr<mi_impactModel> > impactModels_;

  void initializeQP_();

  Eigen::MatrixXd Q_;

  Eigen::MatrixXd C_;

  Eigen::VectorXd cl_, cu_;

  Eigen::VectorXd p_;

  Eigen::VectorXd xl_, xu_;
  Eigen::LSSOL_QP solver_;

  std::vector<std::shared_ptr<mi_equality> > eqConstraints_;

  void solveEqQp_(const Eigen::MatrixXd & Q_,const Eigen::VectorXd & p_, const Eigen::MatrixXd & C_, const Eigen::VectorXd & cu_, Eigen::VectorXd &solution); 


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
