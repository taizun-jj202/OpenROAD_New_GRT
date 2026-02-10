#########
# Setup #
#########

add_library(NEWGR_lib OBJECT
  NewGR.cpp
)

target_include_directories(NEWGR_lib
  PRIVATE
    ..
    include
    ../../include
)

target_link_libraries(NEWGR_lib
  PRIVATE
    utl_lib
    odb
)

add_library(NEWGR_CUGR
  src/CUGR.cpp
  src/Design.cpp
  src/GeoTypes.cpp
  src/GRNet.cpp
  src/GRTree.cpp
  src/GridGraph.cpp
  src/Layers.cpp
  src/Netlist.cpp
  src/MazeRoute.cpp
  src/PatternRoute.cpp
)

target_compile_options(NEWGR_CUGR
  PRIVATE
    -Wno-array-bounds
)

target_include_directories(NEWGR_CUGR
  PUBLIC
    include
  PRIVATE
    # Ugly but necessary to get GRoute.h and avoid a circular dependency
    # with FastRoute.  Remove once this data moves to OpenDB
    ../../include
)

target_link_libraries(NEWGR_CUGR
  PRIVATE
    utl_lib
    stt_lib
    odb
    Boost::boost)
