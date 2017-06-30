cmake_minimum_required(VERSION 3.1)

################################
# axom project
################################
project(axom LANGUAGES C CXX)

if (ENABLE_FORTRAN)
   enable_language(Fortran)
endif()

################################
# BLT
################################
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
    message(FATAL_ERROR "\
The BLT submodule is not present. \
If in git repository run the following two commands:\n \
git submodule init\n \
git submodule update")

endif()

set(ENABLE_COPY_HEADERS ON CACHE PATH "")

include(blt/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)

################################
# Add Thirdparty Builtin Libs
################################
add_subdirectory(thirdparty)

#################################
# Add top-level includes
#################################
add_subdirectory(include)

################################
# Add components
################################
add_subdirectory(components)


################################
# Adds directory for utility code
# that uses or supports axom
################################
add_subdirectory(util)


################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Create header of configuration options
################################s
include(cmake/CMakeConfigureFile.cmake)

################################
# Output a script with all test commands
################################s
if(ENABLE_TESTS)
  set(test_commands_output_file "${PROJECT_BINARY_DIR}/project_test_commands.txt")
  message(STATUS "Outputting the list of test commands to ${test_commands_output_file} ")

  set(test_commands "# List of all test commands in this axom configuration.")
  list(APPEND test_commands "# Format is: test_name  test_command_with_args.")
  list(APPEND test_commands "${BLT_TEST_COMMANDS}")
  string (REPLACE ";" "\n" test_commands_str  "${test_commands}")

  file(WRITE ${test_commands_output_file} "${test_commands_str}" )
endif()

