# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.23

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2022.2.1\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2022.2.1\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\s23049\Documents\GitHub\NAI\CVTest

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\s23049\Documents\GitHub\NAI\CVTest\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/ocvdemo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ocvdemo.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ocvdemo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ocvdemo.dir/flags.make

CMakeFiles/ocvdemo.dir/main.cpp.obj: CMakeFiles/ocvdemo.dir/flags.make
CMakeFiles/ocvdemo.dir/main.cpp.obj: CMakeFiles/ocvdemo.dir/includes_CXX.rsp
CMakeFiles/ocvdemo.dir/main.cpp.obj: ../main.cpp
CMakeFiles/ocvdemo.dir/main.cpp.obj: CMakeFiles/ocvdemo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\s23049\Documents\GitHub\NAI\CVTest\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ocvdemo.dir/main.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.1\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/ocvdemo.dir/main.cpp.obj -MF CMakeFiles\ocvdemo.dir\main.cpp.obj.d -o CMakeFiles\ocvdemo.dir\main.cpp.obj -c C:\Users\s23049\Documents\GitHub\NAI\CVTest\main.cpp

CMakeFiles/ocvdemo.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ocvdemo.dir/main.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.1\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\s23049\Documents\GitHub\NAI\CVTest\main.cpp > CMakeFiles\ocvdemo.dir\main.cpp.i

CMakeFiles/ocvdemo.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ocvdemo.dir/main.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.1\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\s23049\Documents\GitHub\NAI\CVTest\main.cpp -o CMakeFiles\ocvdemo.dir\main.cpp.s

# Object files for target ocvdemo
ocvdemo_OBJECTS = \
"CMakeFiles/ocvdemo.dir/main.cpp.obj"

# External object files for target ocvdemo
ocvdemo_EXTERNAL_OBJECTS =

ocvdemo.exe: CMakeFiles/ocvdemo.dir/main.cpp.obj
ocvdemo.exe: CMakeFiles/ocvdemo.dir/build.make
ocvdemo.exe: C:/Users/s23049/Downloads/opencv-4.x/cmake-build-release/lib/libopencv_highgui460.dll.a
ocvdemo.exe: C:/Users/s23049/Downloads/opencv-4.x/cmake-build-release/lib/libopencv_videoio460.dll.a
ocvdemo.exe: C:/Users/s23049/Downloads/opencv-4.x/cmake-build-release/lib/libopencv_imgcodecs460.dll.a
ocvdemo.exe: C:/Users/s23049/Downloads/opencv-4.x/cmake-build-release/lib/libopencv_imgproc460.dll.a
ocvdemo.exe: C:/Users/s23049/Downloads/opencv-4.x/cmake-build-release/lib/libopencv_core460.dll.a
ocvdemo.exe: CMakeFiles/ocvdemo.dir/linklibs.rsp
ocvdemo.exe: CMakeFiles/ocvdemo.dir/objects1.rsp
ocvdemo.exe: CMakeFiles/ocvdemo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\s23049\Documents\GitHub\NAI\CVTest\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ocvdemo.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\ocvdemo.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ocvdemo.dir/build: ocvdemo.exe
.PHONY : CMakeFiles/ocvdemo.dir/build

CMakeFiles/ocvdemo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\ocvdemo.dir\cmake_clean.cmake
.PHONY : CMakeFiles/ocvdemo.dir/clean

CMakeFiles/ocvdemo.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\s23049\Documents\GitHub\NAI\CVTest C:\Users\s23049\Documents\GitHub\NAI\CVTest C:\Users\s23049\Documents\GitHub\NAI\CVTest\cmake-build-debug C:\Users\s23049\Documents\GitHub\NAI\CVTest\cmake-build-debug C:\Users\s23049\Documents\GitHub\NAI\CVTest\cmake-build-debug\CMakeFiles\ocvdemo.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ocvdemo.dir/depend

