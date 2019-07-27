set(mc_impact_predictor_src
  mi_osd.cpp
  mi_impactPredictor.cpp
  mi_multiImpactPredictor.cpp
  mi_lcp.cpp
  mi_qpEstimator.cpp
)

add_library(mc_impact_predictor SHARED ${mc_impact_predictor_src})
pkg_config_use_dependency(mc_impact_predictor mc_rtc nlopt)
target_link_libraries(mc_impact_predictor ${Boost_LIBRARIES})
target_link_libraries(mc_impact_predictor nlopt)
set_target_properties(mc_impact_predictor PROPERTIES COMPILE_FLAGS "-std=c++0x")

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include/mc_prediction
    FILES_MATCHING PATTERN "*.h*")

install(TARGETS mc_impact_predictor 
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )



