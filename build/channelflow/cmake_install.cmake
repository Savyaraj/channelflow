# Install script for directory: /home/savya/Work/ecps/Scripts/channelflow/channelflow

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libchflow.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libchflow.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libchflow.so"
         RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/savya/Work/ecps/Scripts/channelflow/build/channelflow/libchflow.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libchflow.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libchflow.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libchflow.so"
         OLD_RPATH "/home/savya/Work/ecps/Scripts/channelflow/build/nsolver:/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:"
         NEW_RPATH "/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libchflow.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/channelflow" TYPE FILE FILES
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/bandedtridiag.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/basisfunc.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/cfmpi.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/chebyshev.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/diffops.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/dnsflags.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/dnsalgo.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/nse.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/dns.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/flowfield.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/helmholtz.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/periodicfunc.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/poissonsolver.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/realprofile.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/realprofileng.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/symmetry.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/tausolver.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/turbstats.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/utilfuncs.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/cfdsi.h"
    "/home/savya/Work/ecps/Scripts/channelflow/channelflow/dedalusdsi.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/channelflow" TYPE FILE FILES "/home/savya/Work/ecps/Scripts/channelflow/build/channelflow/config.h")
endif()

