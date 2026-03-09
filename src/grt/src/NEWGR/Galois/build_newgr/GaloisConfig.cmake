# Config file for the Galois package
#
# It exports the following targets:
#   Galois::shmem
#   Galois::dist
#   ...
#   (see GaloisTargets.cmake for all of them)
#
# It defines the following variables for legacy importing:
#   Galois_INCLUDE_DIRS
#   Galois_LIBRARIES
#   Galois_LIBRARY_DIRS
#   Galois_BIN_DIRS
include(CMakeFindDependencyMacro)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was GaloisConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set_and_check(Galois_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set_and_check(Galois_LIBRARY_DIRS "${PACKAGE_PREFIX_DIR}/lib")
set_and_check(Galois_BIN_DIRS "${PACKAGE_PREFIX_DIR}/bin")
set(Galois_LIBRARIES galois_shmem)

find_dependency(Threads REQUIRED)
find_dependency(Boost 1.58.0 REQUIRED COMPONENTS serialization iostreams)
if (OFF)
  find_dependency(MPI REQUIRED)
endif()

get_filename_component(GALOIS_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

if(NOT Galois::shmem)
  include("${GALOIS_CMAKE_DIR}/GaloisTargets.cmake")
endif()
