set(SPROUTE_VENDOR_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/mysproute)

set(GALOIS_ROOT "/home/tsjafri/TAIZUN/Galois" CACHE PATH "Path to Galois source/build tree")

set(GALOIS_INCLUDE_DIRS
    ${GALOIS_ROOT}/libgalois/include
    ${GALOIS_ROOT}/libsupport/include
    ${GALOIS_ROOT}/lonestar/liblonestar/include
)

foreach(dir ${GALOIS_INCLUDE_DIRS})
  if(NOT EXISTS ${dir})
    message(WARNING "Galois include directory not found: ${dir}")
  endif()
endforeach()

file(GLOB SPROUTE_VENDOR_CPP "${SPROUTE_VENDOR_DIR}/src/*.cpp")
file(GLOB SPROUTE_VENDOR_C "${SPROUTE_VENDOR_DIR}/src/*.c")

add_library(SPROUTE_ENGINE STATIC
  SprouteEngine.cpp
  ${SPROUTE_VENDOR_CPP}
  ${SPROUTE_VENDOR_C}
)

target_include_directories(SPROUTE_ENGINE
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${SPROUTE_VENDOR_DIR}/include
    ${GALOIS_INCLUDE_DIRS}
)
