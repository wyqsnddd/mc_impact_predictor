set(mc_impact_predictor_src
  mi_osd.cpp
  mi_lcp.cpp
  mi_qpEstimator.cpp
  mi_jsdEquality.cpp
  mi_invOsdEquality.cpp
  mi_iniEquality.cpp
  mi_impactModel.cpp
  mi_balance.cpp
)

add_library(mc_impact_predictor SHARED ${mc_impact_predictor_src})
target_link_libraries(mc_impact_predictor PUBLIC mc_rtc::mc_rbdyn)
target_link_libraries(mc_impact_predictor PUBLIC eigen-lssol::eigen-lssol)
target_link_libraries(mc_impact_predictor PRIVATE nlopt::nlopt)
target_link_libraries(mc_impact_predictor PUBLIC FIDynamics::FIDynamics)
#target_link_libraries(mc_impact_predictor PRIVATE FIDynamics::FIDynamics)

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



