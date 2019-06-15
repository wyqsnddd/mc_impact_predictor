#include "mi_multiImpactPredictor.h"

mi_multiImpactPredictor::mi_multiImpactPredictor(mc_rbdyn::Robot & robot,
                                                 std::vector<std::string> & impactBodies,
                                                 int maxNumEe,
                                                 double impactDuration,
                                                 double coeFrictionDeduction,
                                                 double coeRes)
: robot_(robot), impactDuration_(impactDuration), coeFrictionDeduction_(coeFrictionDeduction), coeRes_(coeRes)
{

  osdPtr_ = std::make_shared<mi_osd>(getRobot(), true);
  osdPtr_->initializeDataStructure(maxNumEe);
  osdPtr_->resetDataStructure();

  for(auto ii = impactBodies.begin(); ii != impactBodies.end(); ++ii)
  {
    addImpactBody_(*ii);
  }
}

void mi_multiImpactPredictor::run(const std::map<std::string, Eigen::Vector3d> & surfaceNormals)
{
  if(surfaceNormals.size()!=predictorContainer.size() ){
    throw std::runtime_error(
		    std::string("mi_multiImpactPredictor-run: surfaceNormals size(") 
		    + std::to_string(static_cast<int>(surfaceNormals.size()))
		    + std::string(") does not match predictor container size (") 
		    + std::to_string(predictorContainer.size())
		    + std::string(").")
		    );
  }
  osdPtr_->update();

  int dof = getRobot().mb().nrDof();
  jacobianDeltaTau_ = Eigen::MatrixXd::Zero(dof, dof);
  jacobianDeltaAlpha_ = Eigen::MatrixXd::Zero(dof, dof);

  for(auto idx = surfaceNormals.begin(); idx != surfaceNormals.end(); idx++)
  {

    const auto & pi = predictorContainer.find(idx->first);
    pi->second->run(idx->second);

    jacobianDeltaTau_ += pi->second->getJacobianDeltaTau();

    jacobianDeltaAlpha_ += pi->second->getJacobianDeltaAlpha();
  }
  // Update the Jacobians
}

void mi_multiImpactPredictor::addImpactBody_(const std::string & impactBodyName)
{

  const auto & ee = predictorContainer.find(impactBodyName);
  if(ee == predictorContainer.end())
  {
    // create and add
    predictorContainer[impactBodyName] = std::make_shared<mi_impactPredictor>(
        getRobot(), getOsd_(), impactBodyName, true, getImpactDuration(), getCoeFricDe(), getCoeRes());
  }
  else
  {
    throw std::runtime_error("Predictors: Impact body already exists! ");
  }
}

void mi_multiImpactPredictor::setContact(std::vector<std::string> & ees)
{
  for(auto pi = predictorContainer.begin(); pi != predictorContainer.end(); ++pi)
  {
    for(auto index = ees.begin(); index != ees.end(); ++index) pi->second->setContact(*index);
  }
  getOsd_()->setContact(ees);
}

bool mi_multiImpactPredictor::addEndeffectors(const std::vector<std::string> & impactBodies, const std::vector<std::string> & ees)
{
  for (auto impactIdx = impactBodies.begin(); impactIdx!=impactBodies.end(); ++impactIdx){
  // find if the impact exists:
  const auto & pi = predictorContainer.find(*impactIdx);
  if(pi == predictorContainer.end())
  {
    throw std::runtime_error("Predictors-add-ee: Impact body does not exists! ");
  }
  else
  {
    // Initialize the data structure:
    pi->second->initializeDataStructure(static_cast<int>(ees.size()));

    pi->second->resetDataStructure();

    std::cout<<"Initialized impact predictor for impact body: "<<pi->first<<std::endl;
    // Go through all the data
    for(auto index = ees.begin(); index != ees.end(); ++index)
    {
      if(!pi->second->addEndeffector(*index))
      {
        throw std::runtime_error("Impact predictor failed to add endeffector!");
      }
      else
      {
        std::cout << "End-effector: " << *index << " is added to the impact-predictor of: " << pi->first << std::endl;
      }
    } // end of end-effectors
  }

  // add!
 } // end of impact bodies

    return true;
}



const std::shared_ptr<mi_impactPredictor> & mi_multiImpactPredictor::getPredictor(const std::string & impactName)
{
  const auto & pi = predictorContainer.find(impactName);
  if(pi == predictorContainer.end())
  {
    throw std::runtime_error("Predictors-add-ee: Impact body does not exists! ");
  }
  else
    return pi->second;
}
