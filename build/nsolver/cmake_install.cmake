# Install script for directory: /home/savya/Work/ecps/Scripts/channelflow/nsolver

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnsolver.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnsolver.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnsolver.so"
         RPATH "/usr/local/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/savya/Work/ecps/Scripts/channelflow/build/nsolver/libnsolver.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnsolver.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnsolver.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnsolver.so"
         OLD_RPATH "/usr/local/lib:"
         NEW_RPATH "/usr/local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libnsolver.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nsolver" TYPE FILE FILES
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/arnoldi.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/lanczos.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/bicgstab.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/bicgstabl.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/continuation.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/dsi.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/gmres.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/fgmres.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/newtonalgorithm.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/newton.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/nsolver.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/eigenvals.h"
    "/home/savya/Work/ecps/Scripts/channelflow/nsolver/multiShootingDSI.h"
    )
endif()

