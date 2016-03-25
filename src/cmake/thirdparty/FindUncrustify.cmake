################################################################################
# Example usage:
#
# find_package(Uncrustify)
#
# If successful the following variables will be defined
# UNCRUSTIFY_FOUND
# UNCRUSTIFY_EXECUTABLE
################################################################################

find_program(UNCRUSTIFY_EXECUTABLE
             NAMES uncrustify
             DOC "Path to uncrustify executable")

# Handle REQUIRED and QUIET arguments
# this will also set UNCRUSTIFY_FOUND to true if UNCRUSTIFY_EXECUTABLE exists
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Uncrustify
                                  "Failed to locate uncrustify executable"
                                  UNCRUSTIFY_EXECUTABLE)

