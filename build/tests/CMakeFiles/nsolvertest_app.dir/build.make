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
include tests/CMakeFiles/nsolvertest_app.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/nsolvertest_app.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/nsolvertest_app.dir/flags.make

tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o: tests/CMakeFiles/nsolvertest_app.dir/flags.make
tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o: ../tests/nsolvertest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/tests/nsolvertest.cpp

tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/tests/nsolvertest.cpp > CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.i

tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/tests/nsolvertest.cpp -o CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.s

tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.requires:

.PHONY : tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.requires

tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.provides: tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.requires
	$(MAKE) -f tests/CMakeFiles/nsolvertest_app.dir/build.make tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.provides.build
.PHONY : tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.provides

tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.provides.build: tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o


# Object files for target nsolvertest_app
nsolvertest_app_OBJECTS = \
"CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o"

# External object files for target nsolvertest_app
nsolvertest_app_EXTERNAL_OBJECTS =

tests/nsolvertest: tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o
tests/nsolvertest: tests/CMakeFiles/nsolvertest_app.dir/build.make
tests/nsolvertest: channelflow/libchflow.so
tests/nsolvertest: nsolver/libnsolver.so
tests/nsolvertest: /usr/local/lib/libfftw3.so
tests/nsolvertest: tests/CMakeFiles/nsolvertest_app.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable nsolvertest"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nsolvertest_app.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/nsolvertest_app.dir/build: tests/nsolvertest

.PHONY : tests/CMakeFiles/nsolvertest_app.dir/build

tests/CMakeFiles/nsolvertest_app.dir/requires: tests/CMakeFiles/nsolvertest_app.dir/nsolvertest.cpp.o.requires

.PHONY : tests/CMakeFiles/nsolvertest_app.dir/requires

tests/CMakeFiles/nsolvertest_app.dir/clean:
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/nsolvertest_app.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/nsolvertest_app.dir/clean

tests/CMakeFiles/nsolvertest_app.dir/depend:
	cd /home/savya/Work/ecps/Scripts/channelflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/savya/Work/ecps/Scripts/channelflow /home/savya/Work/ecps/Scripts/channelflow/tests /home/savya/Work/ecps/Scripts/channelflow/build /home/savya/Work/ecps/Scripts/channelflow/build/tests /home/savya/Work/ecps/Scripts/channelflow/build/tests/CMakeFiles/nsolvertest_app.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/nsolvertest_app.dir/depend

