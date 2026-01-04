# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

include("openroad")

project(grt)

# Option to switch between FastRoute4.1 and NEWGR
option(USE_NEWGR "Build with the new NEWGR router" ON)
set(USE_NEWGR ON CACHE BOOL "Build with the new NEWGR router" FORCE)

find_package(OpenMP REQUIRED)

if(USE_NEWGR)
  add_subdirectory(src/NEWGR)
else()
  add_subdirectory(src/fastroute)
endif()
add_subdirectory(src/cugr)

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

if(USE_NEWGR)
  target_include_directories(grt_lib
    PRIVATE
      src/NEWGR
      src/NEWGR/include
  )
  target_sources(grt_lib PRIVATE $<TARGET_OBJECTS:NEWGR_lib>)
else()
  target_include_directories(grt_lib
    PRIVATE
      src/fastroute
      src/fastroute/include
  )
endif()

# Conditionally link the router implementation for grt_lib
if(USE_NEWGR)
  message(STATUS "Building grt_lib with NEWGR router")
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
)

if(USE_NEWGR)
  target_sources(grt PRIVATE src/NEWGR/src/FastRouteRenderer.cpp)
else()
  target_sources(grt PRIVATE src/fastroute/src/FastRouteRenderer.cpp)
endif()

target_include_directories(grt
  PUBLIC
    include
  PRIVATE
    src
)

if(USE_NEWGR)
  target_include_directories(grt
    PRIVATE
      src/NEWGR
      src/NEWGR/include
  )
else()
  target_include_directories(grt
    PRIVATE
      src/fastroute
      src/fastroute/include
  )
endif()

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
  if(USE_NEWGR)
    target_include_directories(grt_py
      PRIVATE
        src/NEWGR
        src/NEWGR/include
    )
  else()
    target_include_directories(grt_py
      PRIVATE
        src/fastroute
        src/fastroute/include
    )
  endif()

  # Conditionally link the router implementation for grt_py
  if(USE_NEWGR)
    message(STATUS "Building grt_py with NEWGR router")
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
