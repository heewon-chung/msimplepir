# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.29

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.29.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.29.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build"

# Include any dependencies generated for this target.
include src/CMakeFiles/simplepir.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/simplepir.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/simplepir.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/simplepir.dir/flags.make

src/CMakeFiles/simplepir.dir/database.cpp.o: src/CMakeFiles/simplepir.dir/flags.make
src/CMakeFiles/simplepir.dir/database.cpp.o: /Users/heewonchung/Documents/03.\ Dev/Research/simplepir-cpp/src/database.cpp
src/CMakeFiles/simplepir.dir/database.cpp.o: src/CMakeFiles/simplepir.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/simplepir.dir/database.cpp.o"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/simplepir.dir/database.cpp.o -MF CMakeFiles/simplepir.dir/database.cpp.o.d -o CMakeFiles/simplepir.dir/database.cpp.o -c "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/database.cpp"

src/CMakeFiles/simplepir.dir/database.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simplepir.dir/database.cpp.i"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/database.cpp" > CMakeFiles/simplepir.dir/database.cpp.i

src/CMakeFiles/simplepir.dir/database.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simplepir.dir/database.cpp.s"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/database.cpp" -o CMakeFiles/simplepir.dir/database.cpp.s

src/CMakeFiles/simplepir.dir/parameter.cpp.o: src/CMakeFiles/simplepir.dir/flags.make
src/CMakeFiles/simplepir.dir/parameter.cpp.o: /Users/heewonchung/Documents/03.\ Dev/Research/simplepir-cpp/src/parameter.cpp
src/CMakeFiles/simplepir.dir/parameter.cpp.o: src/CMakeFiles/simplepir.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/simplepir.dir/parameter.cpp.o"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/simplepir.dir/parameter.cpp.o -MF CMakeFiles/simplepir.dir/parameter.cpp.o.d -o CMakeFiles/simplepir.dir/parameter.cpp.o -c "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/parameter.cpp"

src/CMakeFiles/simplepir.dir/parameter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simplepir.dir/parameter.cpp.i"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/parameter.cpp" > CMakeFiles/simplepir.dir/parameter.cpp.i

src/CMakeFiles/simplepir.dir/parameter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simplepir.dir/parameter.cpp.s"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/parameter.cpp" -o CMakeFiles/simplepir.dir/parameter.cpp.s

src/CMakeFiles/simplepir.dir/util.cpp.o: src/CMakeFiles/simplepir.dir/flags.make
src/CMakeFiles/simplepir.dir/util.cpp.o: /Users/heewonchung/Documents/03.\ Dev/Research/simplepir-cpp/src/util.cpp
src/CMakeFiles/simplepir.dir/util.cpp.o: src/CMakeFiles/simplepir.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/simplepir.dir/util.cpp.o"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/simplepir.dir/util.cpp.o -MF CMakeFiles/simplepir.dir/util.cpp.o.d -o CMakeFiles/simplepir.dir/util.cpp.o -c "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/util.cpp"

src/CMakeFiles/simplepir.dir/util.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simplepir.dir/util.cpp.i"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/util.cpp" > CMakeFiles/simplepir.dir/util.cpp.i

src/CMakeFiles/simplepir.dir/util.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simplepir.dir/util.cpp.s"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/util.cpp" -o CMakeFiles/simplepir.dir/util.cpp.s

src/CMakeFiles/simplepir.dir/simplepir.cpp.o: src/CMakeFiles/simplepir.dir/flags.make
src/CMakeFiles/simplepir.dir/simplepir.cpp.o: /Users/heewonchung/Documents/03.\ Dev/Research/simplepir-cpp/src/simplepir.cpp
src/CMakeFiles/simplepir.dir/simplepir.cpp.o: src/CMakeFiles/simplepir.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/simplepir.dir/simplepir.cpp.o"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/simplepir.dir/simplepir.cpp.o -MF CMakeFiles/simplepir.dir/simplepir.cpp.o.d -o CMakeFiles/simplepir.dir/simplepir.cpp.o -c "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/simplepir.cpp"

src/CMakeFiles/simplepir.dir/simplepir.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simplepir.dir/simplepir.cpp.i"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/simplepir.cpp" > CMakeFiles/simplepir.dir/simplepir.cpp.i

src/CMakeFiles/simplepir.dir/simplepir.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simplepir.dir/simplepir.cpp.s"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src/simplepir.cpp" -o CMakeFiles/simplepir.dir/simplepir.cpp.s

# Object files for target simplepir
simplepir_OBJECTS = \
"CMakeFiles/simplepir.dir/database.cpp.o" \
"CMakeFiles/simplepir.dir/parameter.cpp.o" \
"CMakeFiles/simplepir.dir/util.cpp.o" \
"CMakeFiles/simplepir.dir/simplepir.cpp.o"

# External object files for target simplepir
simplepir_EXTERNAL_OBJECTS =

src/libsimplepir.a: src/CMakeFiles/simplepir.dir/database.cpp.o
src/libsimplepir.a: src/CMakeFiles/simplepir.dir/parameter.cpp.o
src/libsimplepir.a: src/CMakeFiles/simplepir.dir/util.cpp.o
src/libsimplepir.a: src/CMakeFiles/simplepir.dir/simplepir.cpp.o
src/libsimplepir.a: src/CMakeFiles/simplepir.dir/build.make
src/libsimplepir.a: src/CMakeFiles/simplepir.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libsimplepir.a"
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && $(CMAKE_COMMAND) -P CMakeFiles/simplepir.dir/cmake_clean_target.cmake
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simplepir.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/simplepir.dir/build: src/libsimplepir.a
.PHONY : src/CMakeFiles/simplepir.dir/build

src/CMakeFiles/simplepir.dir/clean:
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" && $(CMAKE_COMMAND) -P CMakeFiles/simplepir.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/simplepir.dir/clean

src/CMakeFiles/simplepir.dir/depend:
	cd "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp" "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/src" "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build" "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src" "/Users/heewonchung/Documents/03. Dev/Research/simplepir-cpp/build/src/CMakeFiles/simplepir.dir/DependInfo.cmake" "--color=$(COLOR)"
.PHONY : src/CMakeFiles/simplepir.dir/depend

