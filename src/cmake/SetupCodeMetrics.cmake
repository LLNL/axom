#################################################
# Setup code metrics - coverage, profiling, etc
#################################################

################################
# Enable code coverage via gcov
# Note: Only supported for gnu.
################################
if ( (CMAKE_BUILD_TYPE STREQUAL "Debug") AND (ENABLE_CODECOV) )

   if ( (CMAKE_COMPILER_IS_GNUCXX) OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") )
      include(CodeCoverage)
      add_code_coverage_target(coverage make test)
      SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_COVERAGE}" )
      SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_COVERAGE}" )
 	  MESSAGE(STATUS "Code coverage: Detected clang or gnu, compiler flags added.")
   else()
 	   MESSAGE(STATUS "Code coverage: Did not detect clang or gnu, omitting code coverage target.")
   endif()

endif()

find_package(Valgrind)
if (VALGRIND_FOUND)
	set(MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE} CACHE PATH "")
	set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" CACHE PATH "")
endif()

