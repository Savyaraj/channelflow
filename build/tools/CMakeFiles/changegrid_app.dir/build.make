# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/savya/Work/ecps/Scripts/channelflow

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/savya/Work/ecps/Scripts/channelflow/build

# Include any dependencies generated for this target.
include tools/CMakeFiles/changegrid_app.dir/depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/changegrid_app.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/changegrid_app.dir/flags.make

tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o: tools/CMakeFiles/changegrid_app.dir/flags.make
tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o: ../tools/changegrid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/changegrid_app.dir/changegrid.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/tools/changegrid.cpp

tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/changegrid_app.dir/changegrid.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/tools/changegrid.cpp > CMakeFiles/changegrid_app.dir/changegrid.cpp.i

tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/changegrid_app.dir/changegrid.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/tools/changegrid.cpp -o CMakeFiles/changegrid_app.dir/changegrid.cpp.s

tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.requires:

.PHONY : tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.requires

tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.provides: tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.requires
	$(MAKE) -f tools/CMakeFiles/changegrid_app.dir/build.make tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.provides.build
.PHONY : tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.provides

tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.provides.build: tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o


# Object files for target changegrid_app
changegrid_app_OBJECTS = \
"CMakeFiles/changegrid_app.dir/changegrid.cpp.o"

# External object files for target changegrid_app
changegrid_app_EXTERNAL_OBJECTS =

tools/changegrid: tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o
tools/changegrid: tools/CMakeFiles/changegrid_app.dir/build.make
tools/changegrid: channelflow/libchflow.so
tools/changegrid: nsolver/libnsolver.so
tools/changegrid: /usr/local/lib/libfftw3_mpi.so
tools/changegrid: /usr/local/lib/libfftw3.so
tools/changegrid: tools/CMakeFiles/changegrid_app.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable changegrid"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/changegrid_app.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/changegrid_app.dir/build: tools/changegrid

.PHONY : tools/CMakeFiles/changegrid_app.dir/build

tools/CMakeFiles/changegrid_app.dir/requires: tools/CMakeFiles/changegrid_app.dir/changegrid.cpp.o.requires

.PHONY : tools/CMakeFiles/changegrid_app.dir/requires

tools/CMakeFiles/changegrid_app.dir/clean:
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/changegrid_app.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/changegrid_app.dir/clean

tools/CMakeFiles/changegrid_app.dir/depend:
	cd /home/savya/Work/ecps/Scripts/channelflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/savya/Work/ecps/Scripts/channelflow /home/savya/Work/ecps/Scripts/channelflow/tools /home/savya/Work/ecps/Scripts/channelflow/build /home/savya/Work/ecps/Scripts/channelflow/build/tools /home/savya/Work/ecps/Scripts/channelflow/build/tools/CMakeFiles/changegrid_app.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/changegrid_app.dir/depend

