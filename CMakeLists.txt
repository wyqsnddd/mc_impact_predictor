cmake_minimum_required(VERSION 3.8)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

include(cmake/base.cmake)

set(PROJECT_NAME mc_impact_predictor)

set(PROJECT_DESCRIPTION "The impact predictor based on mc_rtc")
set(PROJECT_URL https://gite.lirmm.fr/yuquan/mc_impact_predictor)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)
setup_project()

find_package(mc_rtc REQUIRED)
find_package(nlopt REQUIRED)

include_directories(include)
add_subdirectory(src)

setup_project_finalize()
setup_project_package_finalize()

