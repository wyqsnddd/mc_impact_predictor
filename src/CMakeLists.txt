set(mc_impact_predictor_src
  mi_osd.cpp
  mi_impactPredictor.cpp
  mi_multiImpactPredictor.cpp
  mi_lcp.cpp
  mi_qpEstimator.cpp
  mi_jsdEquality.cpp
  mi_osdEquality.cpp
  mi_invOsdEquality.cpp
  mi_iniEquality.cpp
  mi_impactModel.cpp
)

add_library(mc_impact_predictor SHARED ${mc_impact_predictor_src})
target_link_libraries(mc_impact_predictor mc_rtc::mc_control)
target_link_libraries(mc_impact_predictor eigen-lssol)
pkg_config_use_dependency(mc_impact_predictor nlopt)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include/mc_prediction
    FILES_MATCHING PATTERN "*.h*")

install(TARGETS mc_impact_predictor 
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )



