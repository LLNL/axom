# Find Valgrind.

find_program(VALGRIND_EXECUTABLE
             NAMES valgrind
             DOC "Path to valgrind executable")

# Handle REQUIRED and QUIET arguments
# this will also set VALGRIND_FOUND to true if VALGRIND_EXECUTABLE exists
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Valgrind
                                  "Failed to locate valgrind executable"
                                  VALGRIND_EXECUTABLE)
