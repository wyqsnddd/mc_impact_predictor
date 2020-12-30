set(mc_impact_predictor_src
  mi_osd.cpp
  mi_lcp.cpp
  mi_qpEstimator.cpp
  mi_jsdEquality.cpp
  mi_invOsdEquality.cpp
  mi_iniEquality.cpp
  OneDimImpactModel.cpp
  mi_balance.cpp
  ImpactDynamicsModelInterface.cpp
  McImpactDynamicsModelInterface.cpp
  mi_fid_impulse.cpp
  mi_contactConstraint.cpp
  mi_unilateralContactConstraint.cpp
  mi_velIniEquality.cpp
  mi_velEstimator.cpp
  mi_velBalance.cpp
  mi_velInvOsdEquality.cpp
  mi_velJsdEquality.cpp
  mi_frictionCone.cpp
)

add_library(mc_impact_predictor SHARED ${mc_impact_predictor_src})
target_link_libraries(mc_impact_predictor PUBLIC mc_rtc::mc_rbdyn)
target_link_libraries(mc_impact_predictor PUBLIC eigen-lssol::eigen-lssol)
target_link_libraries(mc_impact_predictor PRIVATE nlopt::nlopt)
target_link_libraries(mc_impact_predictor PUBLIC mc_dynamicStability::mc_dynamicStability)
target_link_libraries(mc_impact_predictor PUBLIC RobotInterface::RobotInterface)
target_link_libraries(mc_impact_predictor PUBLIC FIDynamics::FIDynamics)
target_link_libraries(mc_impact_predictor PRIVATE GSL::gsl)
#target_link_libraries(mc_impact_predictor PRIVATE FIDynamics::FIDynamics)
#target_include_directories(mc_impact_predictor FIDynamics::FIDynamics)
#target_include_directories(mc_impact_predictor PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
#target_include_directories(mc_impact_predictor INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include/mc_prediction
    FILES_MATCHING PATTERN "*.h*")

install(TARGETS mc_impact_predictor 
  EXPORT "${TARGETS_EXPORT_NAME}"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )



