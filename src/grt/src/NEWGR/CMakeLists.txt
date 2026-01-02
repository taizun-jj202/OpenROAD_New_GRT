# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

add_library(NEWGR_lib
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
)

target_link_libraries(NEWGR_lib
  PRIVATE
    est_lib
    utl_lib
    stt_lib
    odb
    Boost::boost
)