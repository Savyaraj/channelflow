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
include tools/CMakeFiles/optphaseshift_app.dir/depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/optphaseshift_app.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/optphaseshift_app.dir/flags.make

tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o: tools/CMakeFiles/optphaseshift_app.dir/flags.make
tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o: ../tools/optphaseshift.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/tools/optphaseshift.cpp

tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/tools/optphaseshift.cpp > CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.i

tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/tools/optphaseshift.cpp -o CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.s

tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.requires:

.PHONY : tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.requires

tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.provides: tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.requires
	$(MAKE) -f tools/CMakeFiles/optphaseshift_app.dir/build.make tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.provides.build
.PHONY : tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.provides

tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.provides.build: tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o


# Object files for target optphaseshift_app
optphaseshift_app_OBJECTS = \
"CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o"

# External object files for target optphaseshift_app
optphaseshift_app_EXTERNAL_OBJECTS =

tools/optphaseshift: tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o
tools/optphaseshift: tools/CMakeFiles/optphaseshift_app.dir/build.make
tools/optphaseshift: channelflow/libchflow.so
tools/optphaseshift: nsolver/libnsolver.so
tools/optphaseshift: /usr/local/lib/libfftw3.so
tools/optphaseshift: tools/CMakeFiles/optphaseshift_app.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable optphaseshift"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/optphaseshift_app.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/optphaseshift_app.dir/build: tools/optphaseshift

.PHONY : tools/CMakeFiles/optphaseshift_app.dir/build

tools/CMakeFiles/optphaseshift_app.dir/requires: tools/CMakeFiles/optphaseshift_app.dir/optphaseshift.cpp.o.requires

.PHONY : tools/CMakeFiles/optphaseshift_app.dir/requires

tools/CMakeFiles/optphaseshift_app.dir/clean:
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/optphaseshift_app.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/optphaseshift_app.dir/clean

tools/CMakeFiles/optphaseshift_app.dir/depend:
	cd /home/savya/Work/ecps/Scripts/channelflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/savya/Work/ecps/Scripts/channelflow /home/savya/Work/ecps/Scripts/channelflow/tools /home/savya/Work/ecps/Scripts/channelflow/build /home/savya/Work/ecps/Scripts/channelflow/build/tools /home/savya/Work/ecps/Scripts/channelflow/build/tools/CMakeFiles/optphaseshift_app.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/optphaseshift_app.dir/depend
