set(SPROUTE_VENDOR_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/mysproute)

set(GALOIS_ROOT "/home/tsjafri/TAIZUN/Galois" CACHE PATH "Path to Galois source/build tree")

set(GALOIS_INCLUDE_DIRS
    ${GALOIS_ROOT}/libgalois/include
    ${GALOIS_ROOT}/libsupport/include
    ${GALOIS_ROOT}/lonestar/liblonestar/include
    ${GALOIS_ROOT}/build/libgalois/include
)

foreach(dir ${GALOIS_INCLUDE_DIRS})
  if(NOT EXISTS ${dir})
    message(WARNING "Galois include directory not found: ${dir}")
  endif()
endforeach()

find_library(GALOIS_SHMEM_LIB
  NAMES galois_shmem
  HINTS
    ${GALOIS_ROOT}/build/libgalois
    ${GALOIS_ROOT}/build
  NO_DEFAULT_PATH)

find_library(GALOIS_SUPPORT_LIB
  NAMES galois_support
  HINTS
    ${GALOIS_ROOT}/build/libsupport
    ${GALOIS_ROOT}/build
  NO_DEFAULT_PATH)

if(NOT GALOIS_SHMEM_LIB)
  message(FATAL_ERROR "Unable to locate galois_shmem. Set GALOIS_ROOT or GALOIS_SHMEM_LIB.")
endif()

if(NOT GALOIS_SUPPORT_LIB)
  message(FATAL_ERROR "Unable to locate galois_support. Set GALOIS_ROOT or GALOIS_SUPPORT_LIB.")
endif()

file(GLOB SPROUTE_VENDOR_CPP "${SPROUTE_VENDOR_DIR}/src/*.cpp")
file(GLOB SPROUTE_VENDOR_C "${SPROUTE_VENDOR_DIR}/src/*.c")

set_source_files_properties(${SPROUTE_VENDOR_C} PROPERTIES LANGUAGE CXX)

add_library(SPROUTE_ENGINE STATIC
  SprouteEngine.cpp
  ${SPROUTE_VENDOR_CPP}
  ${SPROUTE_VENDOR_C}
)

target_include_directories(SPROUTE_ENGINE
  PUBLIC
    ${CMAKE_SOURCE_DIR}/include
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/sta
    ${CMAKE_SOURCE_DIR}/src/sta/include
    ${CMAKE_BINARY_DIR}/src/sta/include
    ${CMAKE_BINARY_DIR}/src/sta/include/sta
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/grt/include
    ${CMAKE_SOURCE_DIR}/src/utl/include
    ${CMAKE_SOURCE_DIR}/src/odb/include
    ${CMAKE_SOURCE_DIR}/src/odb/src/def/def
    ${CMAKE_SOURCE_DIR}/src/odb/src/lef/lef
    ${SPROUTE_VENDOR_DIR}/include
    ${GALOIS_INCLUDE_DIRS}
)

target_link_libraries(SPROUTE_ENGINE
  PRIVATE
    ${GALOIS_SHMEM_LIB}
    ${GALOIS_SUPPORT_LIB}
    Threads::Threads)

target_compile_definitions(SPROUTE_ENGINE
  PRIVATE
    SPROUTE_LUT_DIR="${CMAKE_SOURCE_DIR}/src/stt/src/flt/etc/"
)
