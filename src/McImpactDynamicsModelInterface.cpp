# include "ImpactDynamics/McImpactDynamicsModelInterface.h"

namespace mc_impact
{
McTwoDimModelBridge::McTwoDimModelBridge(
		  const std::shared_ptr<RobotInterface::Robot> robotPtr,
                  const ImpactModelParams & params,
		  const TwoDimModelBridgeParams & bridgeParams)
: TwoDimModelBridge(robotPtr, params, bridgeParams)
{

  std::cout << RoboticsUtils::info << "McTwoDimModelBridge is created." << RoboticsUtils::reset << std::endl;
}

void McTwoDimModelBridge::logImpulseEstimations()
{
  const std::string & bridgeName = getTwoDimModelBridgeParams().name;

  if(getTwoDimModelBridgeParams().gradientApproximation)
  {
    for(auto & result : velCases_)
    {
      logEntries_.emplace_back(bridgeName + "_" + "ImpulseGradientApproximation" + "_" + "ContactVel" + "_"
                               + std::to_string(result.first));
      getHostCtl_()->logger().addLogEntry(logEntries_.back(), [&result]() { return result.second.impulse; });

      logEntries_.emplace_back(bridgeName + "_" + "ImpulseGradientApproximation" + "_" + "comVelJump" + "_"
                               + std::to_string(result.first));
      getHostCtl_()->logger().addLogEntry(logEntries_.back(), [&result]() { return result.second.linearVelJump; });

      logEntries_.emplace_back(bridgeName + "_" + "ImpulseGradientApproximation" + "_" + "angularVelJump" + "_"
                               + std::to_string(result.first));
      getHostCtl_()->logger().addLogEntry(logEntries_.back(), [&result]() { return result.second.anguleVelJump; });
    }
  }

  //logEntries_.emplace_back(bridgeName + "_" + "curveFitting_comVelJump_c0");
  //getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().c0; });

  logEntries_.emplace_back(bridgeName + "_" + "curveFitting_comVelJump_c1");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().c; });

  logEntries_.emplace_back(bridgeName + "_" + "computationTime");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return twoDimModelPtr_->computationTime(); });

  logEntries_.emplace_back(bridgeName + "_" + "rotationAngle");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return rotationAngle_; });

  logEntries_.emplace_back(bridgeName + "_" + "robot_postImpact_linearVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().linearVel; });

  logEntries_.emplace_back(bridgeName + "_" + "robot_linearVelJump");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                      [this]() { return getRobotPostImpactStates().linearVelJump; });

  logEntries_.emplace_back(bridgeName + "_" + "robot_postImpact_angularVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().anguleVel; });

  logEntries_.emplace_back(bridgeName + "_" + "robot_postImpact_angularVelJump");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(),
                                      [this]() { return getRobotPostImpactStates().anguleVelJump; });

  logEntries_.emplace_back(bridgeName + "_" + "robot_postImpact_impulse");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getRobotPostImpactStates().impulse; });

  std::cout << RoboticsUtils::info << "The impulse is: " << getRobotPostImpactStates().impulse << std::endl;

  logEntries_.emplace_back(bridgeName + "_" + "average_angularVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return rAverageAngularVel_; });

  logEntries_.emplace_back(bridgeName + "_" + "average_linearVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return rAverageLinearVel_; });

  logEntries_.emplace_back(bridgeName + "_" + "impactEventSequence");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() {
    switch(twoDimModelPtr_->getSolution().sequence)
    {
      case FIDynamics::EventSequence::not_specified:
        return 0;
      case FIDynamics::EventSequence::C_R_Process:
        return 1;
      case FIDynamics::EventSequence::C_S_R_Process:
        return 2;
      case FIDynamics::EventSequence::S_C_R_Process:
        return 3;
      case FIDynamics::EventSequence::SpecialSolution:
        return 4;
      default:
        RoboticsUtils::throw_runtime_error("impact event sequence is not defined", __FILE__, __LINE__);
    }
  });
  /*
  logEntries_.emplace_back(bridgeName + "_"+ "object_postImpact_linearVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getObjectPostImpactStates().linearVel; });

  logEntries_.emplace_back(bridgeName + "_"+ "object_postImpact_angularVel");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getObjectPostImpactStates().anguleVel; });

  logEntries_.emplace_back(bridgeName + "_"+ "object_postImpact_impulse");
  getHostCtl_()->logger().addLogEntry(logEntries_.back(), [this]() { return getObjectPostImpactStates().impulse; });
*/
}
void McTwoDimModelBridge::removeImpulseEstimations()
{
  assert(getHostCtl_() != nullptr);

  for(auto & name : logEntries_)
  {
    getHostCtl_()->logger().removeLogEntry(name);
  }
}
} // End of namespace mc_impact
