# Install script for directory: /home/savya/Work/ecps/Scripts/channelflow/cfbasics

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/cfbasics" TYPE FILE FILES
    "/home/savya/Work/ecps/Scripts/channelflow/cfbasics/arglist.h"
    "/home/savya/Work/ecps/Scripts/channelflow/cfbasics/brent.h"
    "/home/savya/Work/ecps/Scripts/channelflow/cfbasics/cfarray.h"
    "/home/savya/Work/ecps/Scripts/channelflow/cfbasics/cfbasics.h"
    "/home/savya/Work/ecps/Scripts/channelflow/cfbasics/cfvector.h"
    "/home/savya/Work/ecps/Scripts/channelflow/cfbasics/complexdefs.h"
    "/home/savya/Work/ecps/Scripts/channelflow/cfbasics/mathdefs.h"
    )
endif()

