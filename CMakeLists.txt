# Set the minimum required version of CMake
# FetchContent needs at least 3.11
cmake_minimum_required(VERSION 3.11)
project(SequenceAligner)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find GTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# Make GTest available, but don't build its tests
set(GTEST_DISABLE_TESTS_PROJECT ON) 
FetchContent_MakeAvailable(googletest)

# --- Library Definition ---

# This compiles the NeedlemanWunsch code into a library
add_library(aligner_lib
  src/NeedlemanWunsch.cpp
)

# Tell the library where to find its headers in the 'include/' folder.
# PUBLIC means anyone who links to aligner_lib will also get this include path.
target_include_directories(aligner_lib PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# --- Testing ---

# This enables the 'ctest' command
enable_testing()

# This tells CMake to look in the 'test' folder
# and process the CMakeLists.txt file it finds there.
add_subdirectory(test)

# --- Main Executable ---

# Add an executable target named 'sequence_aligner'
add_executable(sequence_aligner
    src/main.cpp
)

# Link the executable against the library.
# This gives the executable all the code from aligner_lib
# and also gives it the include path from target_include_directories.
target_link_libraries(sequence_aligner PRIVATE aligner_lib)