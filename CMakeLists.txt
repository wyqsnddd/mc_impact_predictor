cmake_minimum_required(VERSION 3.1)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

set(PROJECT_NAME mc_impact_predictor)
set(PROJECT_DESCRIPTION "The impact predictor based on mc_rtc")
set(PROJECT_VERSION 0.1.0)
set(PROJECT_URL https://gite.lirmm.fr/mc-controllers/mc_impact_predictor)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CXX_DISABLE_WERROR 1)
set(CMAKE_CXX_STANDARD 11)
set(DOXYGEN_USE_MATHJAX "YES")

include(cmake/base.cmake)


project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_project_dependency(eigen-lssol REQUIRED)
add_project_dependency(mc_rtc REQUIRED)
add_project_dependency(RobotInterface REQUIRED)
add_project_dependency(mc_dynamicStability REQUIRED)
add_project_dependency(FIDynamics REQUIRED)

find_package(GSL REQUIRED)
find_package(nlopt REQUIRED)

include_directories(include)

add_subdirectory(src)

