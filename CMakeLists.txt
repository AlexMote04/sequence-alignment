cmake_minimum_required(VERSION 3.18)

# 1. Project Setup
# We enable both C++ and CUDA languages
project(NeedlemanWunsch LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)

# Optimization Flags (Host and Device)
if(MSVC)
    add_compile_options(/O2)
else()
    add_compile_options(-O3)
endif()

# 1. Find the CUDA Toolkit
find_package(CUDAToolkit REQUIRED)

# 2. Source Files
# We group the core logic files here
set(NW_SOURCES
    nw_strategy.h
    CudaBasicNW.cu
)

# 3. Main Benchmark Executable
add_executable(nw_benchmark main.cpp ${NW_SOURCES})

# This adds the include paths for <cuda_runtime.h> to main.cpp
target_link_libraries(nw_benchmark PRIVATE CUDA::cudart)

# Tell CMake to use the CUDA compiler for .cu files
set_target_properties(nw_benchmark PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# -----------------------------------------------------------------------------
# 4. Google Test Setup (Automatically downloads GTest)
# -----------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent GTest from overriding our compiler runtime settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# -----------------------------------------------------------------------------
# 5. Test Executable
# -----------------------------------------------------------------------------
add_executable(nw_tests test_nw.cpp ${NW_SOURCES})

# Link GTest
target_link_libraries(nw_tests PRIVATE gtest_main)

# FIX: Link CUDA Runtime to tests as well!
target_link_libraries(nw_tests PRIVATE CUDA::cudart)
set_target_properties(nw_tests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Enable CTest (so you can run "ctest" in terminal)
enable_testing()
add_test(NAME AllTests COMMAND nw_tests)