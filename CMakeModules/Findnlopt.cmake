#
# Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

# Try to find the nlopt library
#
# This can use NLOPT_PREFIX as an hint
#
# Defines the nlopt::nlopt target

set(NLOPT nlopt::nlopt)
if(NOT TARGET ${NLOPT})
  if(NOT DEFINED NLOPT_PREFIX)
    set(NLOPT_PREFIX ${CMAKE_INSTALL_PREFIX})
  endif()

  find_path(NLOPT_INCLUDE_DIR
    NAMES nlopt.h
    HINTS ${NLOPT_PREFIX}
    )
  find_library(NLOPT_LIBRARY
    NAMES nlopt
    PATHS ${NLOPT_PREFIX}
    )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(NLOPT DEFAULT_MSG NLOPT_LIBRARY NLOPT_INCLUDE_DIR)
  mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY)
  if(NLOPT_FOUND)
    add_library(${NLOPT} INTERFACE IMPORTED GLOBAL)
    set_target_properties(${NLOPT} PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${NLOPT_INCLUDE_DIR}
      INTERFACE_LINK_LIBRARIES ${NLOPT_LIBRARY}
      )
  endif()

endif()
