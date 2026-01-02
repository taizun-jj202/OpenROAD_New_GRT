# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

include("openroad")

project(grt)

# Option to switch between FastRoute4.1 and NEWGR
option(USE_NEWGR "Build with the new NEWGR router" OFF)

find_package(OpenMP REQUIRED)

add_subdirectory(src/fastroute)
add_subdirectory(src/cugr)
add_subdirectory(src/NEWGR)

swig_lib(NAME      grt
         NAMESPACE grt
         I_FILE    src/GlobalRouter.i
         SCRIPTS   src/GlobalRouter.tcl
)

add_library(grt_lib
  src/Grid.cpp
  src/Net.cpp
  src/Pin.cpp
  src/RepairAntennas.cpp
  src/Rudy.cpp
  src/GlobalRouter.cpp
  src/SprouteAdapter.cpp
)

target_include_directories(grt_lib
  PUBLIC
    include
  PRIVATE
    src
)

# Conditionally link the router implementation for grt_lib
if(USE_NEWGR)
  message(STATUS "Building grt_lib with NEWGR router")
  target_link_libraries(grt_lib PRIVATE NEWGR_lib)
else()
  message(STATUS "Building grt_lib with default FastRoute4.1 router")
  target_link_libraries(grt_lib PRIVATE FastRoute4.1)
endif()

target_link_libraries(grt_lib
  PRIVATE
    utl_lib
    CUGR
    sproute_tool
    ant_lib
    dpl_lib
    dbSta_lib
    stt_lib
    rsz_lib
    OpenSTA
    Boost::boost
    OpenMP::OpenMP_CXX
)

target_sources(grt
  PRIVATE
    src/MakeGlobalRouter.cpp
    src/heatMap.cpp
    src/heatMapRudy.cpp
    src/GrouteRenderer.cpp
    src/fastroute/src/FastRouteRenderer.cpp
)

target_include_directories(grt
  PUBLIC
    include
  PRIVATE
    src
)

target_link_libraries(grt
    grt_lib
    dbSta
    gui
    stt
    rsz
)

if (Python3_FOUND AND BUILD_PYTHON)
  swig_lib(NAME          grt_py
           NAMESPACE     grt
           LANGUAGE      python
           I_FILE        src/GlobalRouter-py.i
           SWIG_INCLUDES ${PROJECT_SOURCE_DIR}/include/grt
           SCRIPTS       ${CMAKE_CURRENT_BINARY_DIR}/grt_py.py
  )

  target_include_directories(grt_py
    PUBLIC
      include
  )

  # Conditionally link the router implementation for grt_py
  if(USE_NEWGR)
    message(STATUS "Building grt_py with NEWGR router")
    target_link_libraries(grt_py PUBLIC NEWGR_lib)
  else()
    message(STATUS "Building grt_py with default FastRoute4.1 router")
    target_link_libraries(grt_py PUBLIC FastRoute4.1)
  endif()

  target_link_libraries(grt_py
    PUBLIC
      utl_lib
      CUGR
      drt_lib
      ant_lib
      dpl_lib
      dbSta
      gui
      stt
      OpenSTA
      Boost::boost
  )

endif()

messages(
  TARGET grt
)

add_subdirectory(test)