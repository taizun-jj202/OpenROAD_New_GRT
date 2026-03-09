# Install script for directory: /root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/lonestar

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/build_newgr/lonestar/liblonestar/cmake_install.cmake")
  include("/root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/build_newgr/lonestar/tutorial_examples/cmake_install.cmake")
  include("/root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/build_newgr/lonestar/analytics/cmake_install.cmake")
  include("/root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/build_newgr/lonestar/eda/cmake_install.cmake")
  include("/root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/build_newgr/lonestar/mining/cmake_install.cmake")
  include("/root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/build_newgr/lonestar/scientific/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/root/OpenROAD_New_GRT/src/grt/src/NEWGR/Galois/build_newgr/lonestar/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
