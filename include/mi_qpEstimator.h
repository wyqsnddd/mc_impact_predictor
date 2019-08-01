# pragma once 
# include <iomanip>
# include "mi_osd.h" 
# include <math.h>
# include <iostream>
# include <Eigen/StdVector>

# include "mi_osd.h"
# include <mc_rbdyn/Robots.h>
# include <eigen-lssol/LSSOL_QP.h>
# include "mi_jsdEquality.h"
# include "mi_osdEquality.h"
# include "mi_invOsdEquality.h"
# include "mi_iniEquality.h"
# include "mi_impactModel.h"
# include <limits>

class mi_qpEstimator;
struct qpEstimatorParameter{
  double Qweight = 20;
  std::string impactBodyName="r_wrist";
  double impactDuration = 0.005;
  double coeFrictionDeduction = 0.2;
  double coeRes = 0.8;
  int dim = 3;
};

struct endEffector{
 //size_t startingNumber; // Starting number in the vector
 Eigen::VectorXd eeVJump; 
 Eigen::VectorXd estimatedImpulse;
 Eigen::Vector3d estimatedAverageImpulsiveForce;
};

class mi_qpEstimator{
  public: 
  mi_qpEstimator(const mc_rbdyn::Robot & simRobot,
		const std::shared_ptr<mi_osd> & osdPtr,
		const struct qpEstimatorParameter params
		);
  ~mi_qpEstimator(){
 
  }
  void update(const Eigen::Vector3d & surfaceNormal);
  inline const mc_rbdyn::Robot & getSimRobot()
  {
    return simRobot_;
  }
  
  const Eigen::VectorXd  & getPredictedImpulse(const std::string & bodyName);

  
  void addEndeffector(std::string eeName);

 
  inline int getDof() const
  {
   return getOsd_()->getDof(); 
  }

  inline double getQweight() const
  {
   return params_.Qweight; 
  }
  
  
  
  inline const Eigen::VectorXd & getJointVelJump() 
  {
    return jointVelJump_; 
  }

  const endEffector & getEndeffector( const std::string& name);
  void print();
  inline const std::unique_ptr<mi_impactModel> & getImpactModel() const
  {
    return impactModelPtr_; 
  }
  private: 
  const mc_rbdyn::Robot & simRobot_;
  const std::shared_ptr<mi_osd> & osdPtr_;
  inline const std::shared_ptr<mi_osd> & getOsd_() const
  {
    return osdPtr_;
  }
  endEffector & getEndeffector_( const std::string& name);
  qpEstimatorParameter params_;
  std::unique_ptr<mi_impactModel>  impactModelPtr_;

  void initializeQP_();


  Eigen::MatrixXd Q_;

  Eigen::MatrixXd C_;

  Eigen::VectorXd cl_, cu_;

  Eigen::VectorXd p_;

  Eigen::VectorXd xl_, xu_;
  Eigen::LSSOL_QP solver_;

  std::vector<std::shared_ptr<mi_equality> > eqConstraints_;



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
  Eigen::VectorXd jointVelJump_; 
  
};
