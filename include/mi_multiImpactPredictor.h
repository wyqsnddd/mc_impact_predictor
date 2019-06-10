#pragma once

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/Robots.h>

#include "mi_impactPredictor.h"
#include <map>

class mi_multiImpactPredictor
/** \brief This is a wrapper of a group of mi_impactPredictor.  
 */
{
public:
/*!
      \param robot reference to the robot model used by the QP
      \param impactBodies that expect an impact
      \param maxNumEe the maximal number of end-effectors of each predictor. We use it to allocate memory for the predictors.  
      \param coeFrictionDeduction tangential velocity deduction coeffecient at the impact.
      \param coeRes coefficient of restitution
      */
  mi_multiImpactPredictor(mc_rbdyn::Robot & robot,
                          std::vector<std::string> & impactBodies,
                          int maxNumEe,
                          double impactDuration,
                          double coeFrictionDeduction,
                          double coeRes = 0.8);

  ~mi_multiImpactPredictor() {}

  void run(const std::map<std::string, Eigen::Vector3d> & surfaceNormals);

  inline void resetDataStructure()
  {
    osdPtr_->resetDataStructure();
    for(auto ii = predictorContainer.begin(); ii != predictorContainer.end(); ii++) ii->second->resetDataStructure();
  }
  inline void initializeDataStructure(int numEE)
  {

    osdPtr_->initializeDataStructure(numEE);
    for(auto ii = predictorContainer.begin(); ii != predictorContainer.end(); ii++)
      ii->second->initializeDataStructure(numEE);
  }

  inline mc_rbdyn::Robot & getRobot()
  {
    return robot_;
  }
  inline const double & getImpactDuration() const
  {
    return impactDuration_;
  }
  inline const double & getCoeRes() const
  {
    return coeRes_;
  }
  inline const double & getCoeFricDe() const
  {
    return coeFrictionDeduction_;
  }
  bool addEndeffectors(const std::vector<std::string> & impactBodies, const std::vector<std::string> & ees);

  /*!
     \param ees includes the bodies with established contacts.
     We set them in contact for all the predictors.
     */
  void setContact(std::vector<std::string> & ees);
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
    return jacobianDeltaTau_;
  }
  inline const std::shared_ptr<mi_osd> & getOsd_() const
  {
    return osdPtr_;
  }
  /*!
   * \param impactName get the impact predictor corresponds to the impactName
   * \brief get a reference to the shared pointer of the impact predictor.
   */
  const std::shared_ptr<mi_impactPredictor> & getPredictor(const std::string & impactName);

protected:
  mc_rbdyn::Robot & robot_;
  std::map<std::string, std::shared_ptr<mi_impactPredictor>> predictorContainer;
  double impactDuration_;
  double coeFrictionDeduction_;
  double coeRes_;

  std::shared_ptr<mi_osd> osdPtr_;
  Eigen::MatrixXd jacobianDeltaAlpha_;
  Eigen::MatrixXd jacobianDeltaTau_;

  void addImpactBody_(const std::string & impactBodyName);
};
