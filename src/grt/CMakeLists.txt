# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

include("openroad")

project(grt)

# Option to build the experimental NEWGR router alongside the legacy
# FastRoute4.1 implementation.
option(USE_NEWGR "Build with the new NEWGR router" ON)

find_package(OpenMP REQUIRED)

add_subdirectory(src/fastroute)
# Always build the canonical CUGR implementation.
add_subdirectory(src/cugr)
# Optionally build NEWGR (and its own CUGR-derived engine) alongside it.
if(USE_NEWGR)
  add_subdirectory(src/NEWGR)
endif()

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

target_include_directories(grt_lib
  PRIVATE
    src/fastroute
    src/fastroute/include
)

if(USE_NEWGR)
  target_include_directories(grt_lib
    PRIVATE
      src/NEWGR
      src/NEWGR/include
  )
  target_sources(grt_lib PRIVATE $<TARGET_OBJECTS:NEWGR_lib>)
  target_compile_definitions(grt_lib PUBLIC GRT_ENABLE_NEWGR=1)
else()
  target_compile_definitions(grt_lib PUBLIC GRT_ENABLE_NEWGR=0)
endif()

set(GRT_EXTRA_LIBS "")
if(USE_NEWGR)
  list(APPEND GRT_EXTRA_LIBS NEWGR_CUGR)
endif()

target_link_libraries(grt_lib
    utl_lib
    FastRoute4.1
    CUGR
    ${GRT_EXTRA_LIBS}
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

target_include_directories(grt
  PUBLIC
    include
  PRIVATE
    src
    src/fastroute
    src/fastroute/include
)

target_sources(grt PRIVATE src/fastroute/src/FastRouteRenderer.cpp)
if(USE_NEWGR)
  target_include_directories(grt
    PRIVATE
      src/NEWGR
      src/NEWGR/include
  )
endif()

target_link_libraries(grt
    grt_lib
    dbSta
    gui
    stt
    rsz
)
if(USE_NEWGR)
  target_compile_definitions(grt PRIVATE GRT_ENABLE_NEWGR=1)
else()
  target_compile_definitions(grt PRIVATE GRT_ENABLE_NEWGR=0)
endif()

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
    PRIVATE
      src/fastroute
      src/fastroute/include
  )

  target_link_libraries(grt_py
    PUBLIC
      utl_lib
      FastRoute4.1
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

  if(USE_NEWGR)
    target_include_directories(grt_py
      PRIVATE
        src/NEWGR
        src/NEWGR/include
    )
    target_compile_definitions(grt_py PRIVATE GRT_ENABLE_NEWGR=1)
  else()
    target_compile_definitions(grt_py PRIVATE GRT_ENABLE_NEWGR=0)
  endif()

endif()

messages(
  TARGET grt
)

add_subdirectory(test)
