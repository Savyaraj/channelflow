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
include nsolver/CMakeFiles/nsolver.dir/depend.make

# Include the progress variables for this target.
include nsolver/CMakeFiles/nsolver.dir/progress.make

# Include the compile flags for this target's objects.
include nsolver/CMakeFiles/nsolver.dir/flags.make

nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o: ../nsolver/lanczos.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/lanczos.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/lanczos.cpp

nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/lanczos.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/lanczos.cpp > CMakeFiles/nsolver.dir/lanczos.cpp.i

nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/lanczos.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/lanczos.cpp -o CMakeFiles/nsolver.dir/lanczos.cpp.s

nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o


nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o: ../nsolver/arnoldi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/arnoldi.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/arnoldi.cpp

nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/arnoldi.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/arnoldi.cpp > CMakeFiles/nsolver.dir/arnoldi.cpp.i

nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/arnoldi.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/arnoldi.cpp -o CMakeFiles/nsolver.dir/arnoldi.cpp.s

nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o


nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o: ../nsolver/bicgstab.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/bicgstab.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/bicgstab.cpp

nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/bicgstab.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/bicgstab.cpp > CMakeFiles/nsolver.dir/bicgstab.cpp.i

nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/bicgstab.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/bicgstab.cpp -o CMakeFiles/nsolver.dir/bicgstab.cpp.s

nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o


nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o: ../nsolver/continuation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/continuation.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/continuation.cpp

nsolver/CMakeFiles/nsolver.dir/continuation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/continuation.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/continuation.cpp > CMakeFiles/nsolver.dir/continuation.cpp.i

nsolver/CMakeFiles/nsolver.dir/continuation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/continuation.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/continuation.cpp -o CMakeFiles/nsolver.dir/continuation.cpp.s

nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o


nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o: ../nsolver/dsi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/dsi.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/dsi.cpp

nsolver/CMakeFiles/nsolver.dir/dsi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/dsi.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/dsi.cpp > CMakeFiles/nsolver.dir/dsi.cpp.i

nsolver/CMakeFiles/nsolver.dir/dsi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/dsi.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/dsi.cpp -o CMakeFiles/nsolver.dir/dsi.cpp.s

nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o


nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o: ../nsolver/gmres.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/gmres.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/gmres.cpp

nsolver/CMakeFiles/nsolver.dir/gmres.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/gmres.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/gmres.cpp > CMakeFiles/nsolver.dir/gmres.cpp.i

nsolver/CMakeFiles/nsolver.dir/gmres.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/gmres.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/gmres.cpp -o CMakeFiles/nsolver.dir/gmres.cpp.s

nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o


nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o: ../nsolver/fgmres.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/fgmres.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/fgmres.cpp

nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/fgmres.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/fgmres.cpp > CMakeFiles/nsolver.dir/fgmres.cpp.i

nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/fgmres.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/fgmres.cpp -o CMakeFiles/nsolver.dir/fgmres.cpp.s

nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o


nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o: ../nsolver/newtonalgorithm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/newtonalgorithm.cpp

nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/newtonalgorithm.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/newtonalgorithm.cpp > CMakeFiles/nsolver.dir/newtonalgorithm.cpp.i

nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/newtonalgorithm.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/newtonalgorithm.cpp -o CMakeFiles/nsolver.dir/newtonalgorithm.cpp.s

nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o


nsolver/CMakeFiles/nsolver.dir/newton.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/newton.cpp.o: ../nsolver/newton.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object nsolver/CMakeFiles/nsolver.dir/newton.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/newton.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/newton.cpp

nsolver/CMakeFiles/nsolver.dir/newton.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/newton.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/newton.cpp > CMakeFiles/nsolver.dir/newton.cpp.i

nsolver/CMakeFiles/nsolver.dir/newton.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/newton.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/newton.cpp -o CMakeFiles/nsolver.dir/newton.cpp.s

nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/newton.cpp.o


nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o: ../nsolver/eigenvals.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/eigenvals.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/eigenvals.cpp

nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/eigenvals.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/eigenvals.cpp > CMakeFiles/nsolver.dir/eigenvals.cpp.i

nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/eigenvals.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/eigenvals.cpp -o CMakeFiles/nsolver.dir/eigenvals.cpp.s

nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o


nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o: nsolver/CMakeFiles/nsolver.dir/flags.make
nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o: ../nsolver/multiShootingDSI.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o -c /home/savya/Work/ecps/Scripts/channelflow/nsolver/multiShootingDSI.cpp

nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nsolver.dir/multiShootingDSI.cpp.i"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/savya/Work/ecps/Scripts/channelflow/nsolver/multiShootingDSI.cpp > CMakeFiles/nsolver.dir/multiShootingDSI.cpp.i

nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nsolver.dir/multiShootingDSI.cpp.s"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && /usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/savya/Work/ecps/Scripts/channelflow/nsolver/multiShootingDSI.cpp -o CMakeFiles/nsolver.dir/multiShootingDSI.cpp.s

nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.requires:

.PHONY : nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.requires

nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.provides: nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.requires
	$(MAKE) -f nsolver/CMakeFiles/nsolver.dir/build.make nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.provides.build
.PHONY : nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.provides

nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.provides.build: nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o


# Object files for target nsolver
nsolver_OBJECTS = \
"CMakeFiles/nsolver.dir/lanczos.cpp.o" \
"CMakeFiles/nsolver.dir/arnoldi.cpp.o" \
"CMakeFiles/nsolver.dir/bicgstab.cpp.o" \
"CMakeFiles/nsolver.dir/continuation.cpp.o" \
"CMakeFiles/nsolver.dir/dsi.cpp.o" \
"CMakeFiles/nsolver.dir/gmres.cpp.o" \
"CMakeFiles/nsolver.dir/fgmres.cpp.o" \
"CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o" \
"CMakeFiles/nsolver.dir/newton.cpp.o" \
"CMakeFiles/nsolver.dir/eigenvals.cpp.o" \
"CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o"

# External object files for target nsolver
nsolver_EXTERNAL_OBJECTS =

nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/newton.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/build.make
nsolver/libnsolver.so: /usr/local/lib/libfftw3.so
nsolver/libnsolver.so: nsolver/CMakeFiles/nsolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/savya/Work/ecps/Scripts/channelflow/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX shared library libnsolver.so"
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nsolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
nsolver/CMakeFiles/nsolver.dir/build: nsolver/libnsolver.so

.PHONY : nsolver/CMakeFiles/nsolver.dir/build

nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/lanczos.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/arnoldi.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/bicgstab.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/continuation.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/dsi.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/gmres.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/fgmres.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/newtonalgorithm.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/newton.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/eigenvals.cpp.o.requires
nsolver/CMakeFiles/nsolver.dir/requires: nsolver/CMakeFiles/nsolver.dir/multiShootingDSI.cpp.o.requires

.PHONY : nsolver/CMakeFiles/nsolver.dir/requires

nsolver/CMakeFiles/nsolver.dir/clean:
	cd /home/savya/Work/ecps/Scripts/channelflow/build/nsolver && $(CMAKE_COMMAND) -P CMakeFiles/nsolver.dir/cmake_clean.cmake
.PHONY : nsolver/CMakeFiles/nsolver.dir/clean

nsolver/CMakeFiles/nsolver.dir/depend:
	cd /home/savya/Work/ecps/Scripts/channelflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/savya/Work/ecps/Scripts/channelflow /home/savya/Work/ecps/Scripts/channelflow/nsolver /home/savya/Work/ecps/Scripts/channelflow/build /home/savya/Work/ecps/Scripts/channelflow/build/nsolver /home/savya/Work/ecps/Scripts/channelflow/build/nsolver/CMakeFiles/nsolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : nsolver/CMakeFiles/nsolver.dir/depend

