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
include programs/CMakeFiles/simulateflow_app.dir/depend.make

# Include the progress variables for this target.
include programs/CMakeFiles/simulateflow_app.dir/progress.make

# Include the compile flags for this target's objects.
include programs/CMakeFiles/simulateflow_app.dir/flags.make

programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o: programs/CMakeFiles/simulateflow_app.dir/flags.make
programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o: ../programs/simulateflow.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/programs && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/programs/simulateflow.cpp

programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simulateflow_app.dir/simulateflow.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/programs && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/programs/simulateflow.cpp > CMakeFiles/simulateflow_app.dir/simulateflow.cpp.i

programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simulateflow_app.dir/simulateflow.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/programs && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/programs/simulateflow.cpp -o CMakeFiles/simulateflow_app.dir/simulateflow.cpp.s

programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.requires:

.PHONY : programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.requires

programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.provides: programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.requires
	$(MAKE) -f programs/CMakeFiles/simulateflow_app.dir/build.make programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.provides.build
.PHONY : programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.provides

programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.provides.build: programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o


# Object files for target simulateflow_app
simulateflow_app_OBJECTS = \
"CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o"

# External object files for target simulateflow_app
simulateflow_app_EXTERNAL_OBJECTS =

programs/simulateflow: programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o
programs/simulateflow: programs/CMakeFiles/simulateflow_app.dir/build.make
programs/simulateflow: channelflow/libchflow.so
programs/simulateflow: nsolver/libnsolver.so
programs/simulateflow: /usr/local/lib/libfftw3.so
programs/simulateflow: programs/CMakeFiles/simulateflow_app.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable simulateflow"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/programs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simulateflow_app.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
programs/CMakeFiles/simulateflow_app.dir/build: programs/simulateflow

.PHONY : programs/CMakeFiles/simulateflow_app.dir/build

programs/CMakeFiles/simulateflow_app.dir/requires: programs/CMakeFiles/simulateflow_app.dir/simulateflow.cpp.o.requires

.PHONY : programs/CMakeFiles/simulateflow_app.dir/requires

programs/CMakeFiles/simulateflow_app.dir/clean:
	cd /home/savya/Work/ecps/Scripts/channelflow/build/programs && $(CMAKE_COMMAND) -P CMakeFiles/simulateflow_app.dir/cmake_clean.cmake
.PHONY : programs/CMakeFiles/simulateflow_app.dir/clean

programs/CMakeFiles/simulateflow_app.dir/depend:
	cd /home/savya/Work/ecps/Scripts/channelflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/savya/Work/ecps/Scripts/channelflow /home/savya/Work/ecps/Scripts/channelflow/programs /home/savya/Work/ecps/Scripts/channelflow/build /home/savya/Work/ecps/Scripts/channelflow/build/programs /home/savya/Work/ecps/Scripts/channelflow/build/programs/CMakeFiles/simulateflow_app.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : programs/CMakeFiles/simulateflow_app.dir/depend

