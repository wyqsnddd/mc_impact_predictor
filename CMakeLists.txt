cmake_minimum_required(VERSION 2.8)

include(cmake/base.cmake)
include(cmake/boost.cmake)

set(PROJECT_NAME mc_impact_predictor)

set(PROJECT_DESCRIPTION "The impact predictor based on mc_rtc")
set(PROJECT_URL https://gite.lirmm.fr/yuquan/mc_impact_predictor)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CXX_DISABLE_WERROR 1)
setup_project()

# Enable C++11
if(NOT WIN32)
  add_definitions(-std=c++0x)
endif()

# Use jrl-cmakemoduels to get mc_rtc depency
set(PKG_CONFIG_ADDITIONAL_VARIABLES "${PKG_CONFIG_ADDITIONAL_VARIABLES} libdir")
#add_required_dependency(mc_rtc)
find_package(DART 6.8.0 REQUIRED COMPONENTS utils-urdf gui-osg collision-bullet CONFIG)
find_package(Bullet REQUIRED)
FIND_PACKAGE(Boost  COMPONENTS program_options REQUIRED )


# add_subdirectory(src/example_one)

include_directories(include)
add_subdirectory(src)

setup_project_finalize()

