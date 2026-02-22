set(NEWGR_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(NEWGR_VENDOR_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/mysproute)

# Pin Galois to the in-tree copy for reproducibility.
set(NEWGR_GALOIS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/Galois")
set(NEWGR_GALOIS_BUILD_DIR "${NEWGR_GALOIS_ROOT}/build_newgr")

set(NEWGR_GALOIS_INCLUDE_DIRS
    ${NEWGR_GALOIS_ROOT}/libgalois/include
    ${NEWGR_GALOIS_ROOT}/libsupport/include
    ${NEWGR_GALOIS_ROOT}/lonestar/liblonestar/include
    ${NEWGR_GALOIS_ROOT}/build/libgalois/include
)

foreach(dir ${NEWGR_GALOIS_INCLUDE_DIRS})
  if(NOT EXISTS ${dir})
    message(WARNING "NEWGR Galois include directory not found: ${dir}")
  endif()
endforeach()

find_package(Threads REQUIRED)

find_library(NEWGR_GALOIS_SHMEM_LIB
  NAMES galois_shmem
  HINTS
    ${NEWGR_GALOIS_BUILD_DIR}/libgalois
    ${NEWGR_GALOIS_BUILD_DIR}
    ${NEWGR_GALOIS_ROOT}/build/libgalois
    ${NEWGR_GALOIS_ROOT}/build
  NO_DEFAULT_PATH)

find_library(NEWGR_GALOIS_SUPPORT_LIB
  NAMES galois_support
  HINTS
    ${NEWGR_GALOIS_BUILD_DIR}/libsupport
    ${NEWGR_GALOIS_BUILD_DIR}
    ${NEWGR_GALOIS_ROOT}/build/libsupport
    ${NEWGR_GALOIS_ROOT}/build
  NO_DEFAULT_PATH)

if(NOT NEWGR_GALOIS_SHMEM_LIB)
  message(FATAL_ERROR "Unable to locate NEWGR galois_shmem. Set NEWGR_GALOIS_ROOT or NEWGR_GALOIS_SHMEM_LIB.")
endif()

if(NOT NEWGR_GALOIS_SUPPORT_LIB)
  message(FATAL_ERROR "Unable to locate NEWGR galois_support. Set NEWGR_GALOIS_ROOT or NEWGR_GALOIS_SUPPORT_LIB.")
endif()

file(GLOB NEWGR_VENDOR_CPP "${NEWGR_VENDOR_DIR}/src/*.cpp")

add_library(newgr_tool SHARED
  ${NEWGR_ROOT}/NewgrEngine.cpp
  ${NEWGR_VENDOR_CPP}
)

target_link_options(newgr_tool
  PRIVATE
    -Wl,-Bsymbolic
)

target_include_directories(newgr_tool
  PRIVATE
    ${NEWGR_ROOT}
    ${NEWGR_VENDOR_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/grt/include
    ${CMAKE_SOURCE_DIR}/src/utl/include
    ${CMAKE_SOURCE_DIR}/src/odb/include
    ${CMAKE_SOURCE_DIR}/src/odb/src/def/def
    ${CMAKE_SOURCE_DIR}/src/odb/src/lef/lef
    ${CMAKE_SOURCE_DIR}/src/sta
    ${CMAKE_SOURCE_DIR}/src/sta/include
    ${CMAKE_BINARY_DIR}/src/sta/include
    ${CMAKE_BINARY_DIR}/src/sta/include/sta
    ${NEWGR_GALOIS_INCLUDE_DIRS}
)

target_link_libraries(newgr_tool
  PRIVATE
    ${NEWGR_GALOIS_SHMEM_LIB}
    ${NEWGR_GALOIS_SUPPORT_LIB}
    Threads::Threads
)

target_compile_definitions(newgr_tool
  PRIVATE
    SPROUTE_LUT_DIR="${CMAKE_SOURCE_DIR}/src/stt/src/flt/etc/"
)
