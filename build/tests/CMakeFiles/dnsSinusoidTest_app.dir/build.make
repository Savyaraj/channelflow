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
include tests/CMakeFiles/dnsSinusoidTest_app.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/dnsSinusoidTest_app.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/dnsSinusoidTest_app.dir/flags.make

tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o: tests/CMakeFiles/dnsSinusoidTest_app.dir/flags.make
tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o: ../tests/dnsSinusoidTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/tests/dnsSinusoidTest.cpp

tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/tests/dnsSinusoidTest.cpp > CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.i

tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/tests/dnsSinusoidTest.cpp -o CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.s

tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.requires:

.PHONY : tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.requires

tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.provides: tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.requires
	$(MAKE) -f tests/CMakeFiles/dnsSinusoidTest_app.dir/build.make tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.provides.build
.PHONY : tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.provides

tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.provides.build: tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o


# Object files for target dnsSinusoidTest_app
dnsSinusoidTest_app_OBJECTS = \
"CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o"

# External object files for target dnsSinusoidTest_app
dnsSinusoidTest_app_EXTERNAL_OBJECTS =

tests/dnsSinusoidTest: tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o
tests/dnsSinusoidTest: tests/CMakeFiles/dnsSinusoidTest_app.dir/build.make
tests/dnsSinusoidTest: channelflow/libchflow.so
tests/dnsSinusoidTest: nsolver/libnsolver.so
tests/dnsSinusoidTest: /usr/local/lib/libfftw3_mpi.so
tests/dnsSinusoidTest: /usr/local/lib/libfftw3.so
tests/dnsSinusoidTest: tests/CMakeFiles/dnsSinusoidTest_app.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dnsSinusoidTest"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dnsSinusoidTest_app.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/dnsSinusoidTest_app.dir/build: tests/dnsSinusoidTest

.PHONY : tests/CMakeFiles/dnsSinusoidTest_app.dir/build

tests/CMakeFiles/dnsSinusoidTest_app.dir/requires: tests/CMakeFiles/dnsSinusoidTest_app.dir/dnsSinusoidTest.cpp.o.requires

.PHONY : tests/CMakeFiles/dnsSinusoidTest_app.dir/requires

tests/CMakeFiles/dnsSinusoidTest_app.dir/clean:
	cd /home/savya/Work/ecps/Scripts/channelflow/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/dnsSinusoidTest_app.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/dnsSinusoidTest_app.dir/clean

tests/CMakeFiles/dnsSinusoidTest_app.dir/depend:
	cd /home/savya/Work/ecps/Scripts/channelflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/savya/Work/ecps/Scripts/channelflow /home/savya/Work/ecps/Scripts/channelflow/tests /home/savya/Work/ecps/Scripts/channelflow/build /home/savya/Work/ecps/Scripts/channelflow/build/tests /home/savya/Work/ecps/Scripts/channelflow/build/tests/CMakeFiles/dnsSinusoidTest_app.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/dnsSinusoidTest_app.dir/depend

