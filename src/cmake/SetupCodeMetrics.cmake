#################################################
# Setup code metrics - coverage, profiling, etc
#################################################

########################################
# Enable code coverage via gcov
# Note: Only supported for gnu or clang.
########################################
if (ENABLE_CODECOV)
   include(CodeCoverage)
endif()

find_package(Valgrind)
if (VALGRIND_FOUND)
	set(MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE} CACHE PATH "")
	set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" CACHE PATH "")
endif()

