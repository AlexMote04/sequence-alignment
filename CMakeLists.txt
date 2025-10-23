# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)
project(SequenceAligner)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add an executable target named 'sequence_aligner'
# and list all the .cpp files needed to build it.
add_executable(sequence_aligner
    src/main.cpp
    src/NeedlemanWunsch.cpp
)

# Tell the compiler where to find the header files
target_include_directories(sequence_aligner PUBLIC include)