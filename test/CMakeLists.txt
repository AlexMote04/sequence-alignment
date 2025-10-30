# Define the Test Executable
add_executable(run_tests
  test_main.cpp
)

# Link Libraries
target_link_libraries(run_tests PRIVATE
  aligner_lib  
  gtest_main   # This provides the main() function for you
  gmock        # This is needed for UnorderedElementsAre
)

# Register the Tests with CTest
# This finds all the TEST_F/TEST P macros
# and adds them to the 'ctest' runner.
include(GoogleTest)
gtest_discover_tests(run_tests)