# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

# Compute the OpenROAD repository root relative to this CMake file so the
# include paths remain valid even if OpenROAD is built as a sub-project.
get_filename_component(NEWGR_OPENROAD_ROOT
                       "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)

add_library(NEWGR_lib OBJECT
  src/FastRoute.cpp
  src/RSMT.cpp
  src/RipUp.cpp
  src/maze.cpp
  src/maze3D.cpp
  src/route.cpp
  src/utility.cpp
  src/graph2d.cpp
)

target_compile_options(NEWGR_lib
  PRIVATE
    -Wno-array-bounds
)

target_include_directories(NEWGR_lib
  PUBLIC
    include
  PRIVATE
    # Ugly but necessary to get GRoute.h and avoid a circular dependency
    # with FastRoute.  Remove once this data moves to OpenDB
    ../../include
    ${NEWGR_OPENROAD_ROOT}/src
    ${NEWGR_OPENROAD_ROOT}/src/utl/include
    ${NEWGR_OPENROAD_ROOT}/src/odb/include
    ${NEWGR_OPENROAD_ROOT}/src/dbSta/include
    ${NEWGR_OPENROAD_ROOT}/src/stt/include
    ${NEWGR_OPENROAD_ROOT}/src/sta/include
    ${CMAKE_BINARY_DIR}/src/sta/include
    ${CMAKE_BINARY_DIR}/src/sta/include/sta
    ..
)
