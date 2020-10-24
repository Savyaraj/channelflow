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
include tools/CMakeFiles/findsymmetries_app.dir/depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/findsymmetries_app.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/findsymmetries_app.dir/flags.make

tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o: tools/CMakeFiles/findsymmetries_app.dir/flags.make
tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o: ../tools/findsymmetries.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/tools/findsymmetries.cpp

tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/tools/findsymmetries.cpp > CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.i

tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/tools/findsymmetries.cpp -o CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.s

tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.requires:

.PHONY : tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.requires

tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.provides: tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.requires
	$(MAKE) -f tools/CMakeFiles/findsymmetries_app.dir/build.make tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.provides.build
.PHONY : tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.provides

tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.provides.build: tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o


# Object files for target findsymmetries_app
findsymmetries_app_OBJECTS = \
"CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o"

# External object files for target findsymmetries_app
findsymmetries_app_EXTERNAL_OBJECTS =

tools/findsymmetries: tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o
tools/findsymmetries: tools/CMakeFiles/findsymmetries_app.dir/build.make
tools/findsymmetries: channelflow/libchflow.so
tools/findsymmetries: nsolver/libnsolver.so
tools/findsymmetries: /usr/local/lib/libfftw3.so
tools/findsymmetries: tools/CMakeFiles/findsymmetries_app.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable findsymmetries"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/findsymmetries_app.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/findsymmetries_app.dir/build: tools/findsymmetries

.PHONY : tools/CMakeFiles/findsymmetries_app.dir/build

tools/CMakeFiles/findsymmetries_app.dir/requires: tools/CMakeFiles/findsymmetries_app.dir/findsymmetries.cpp.o.requires

.PHONY : tools/CMakeFiles/findsymmetries_app.dir/requires

tools/CMakeFiles/findsymmetries_app.dir/clean:
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/findsymmetries_app.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/findsymmetries_app.dir/clean

tools/CMakeFiles/findsymmetries_app.dir/depend:
	cd /home/savya/Work/ecps/Scripts/channelflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/savya/Work/ecps/Scripts/channelflow /home/savya/Work/ecps/Scripts/channelflow/tools /home/savya/Work/ecps/Scripts/channelflow/build /home/savya/Work/ecps/Scripts/channelflow/build/tools /home/savya/Work/ecps/Scripts/channelflow/build/tools/CMakeFiles/findsymmetries_app.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/findsymmetries_app.dir/depend

