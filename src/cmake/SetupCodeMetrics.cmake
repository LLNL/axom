#################################################
# Setup code metrics - coverage, profiling, etc
#################################################

################################
# Enable code coverage via gcov
# Note: Only supported for gnu.
################################
include(CodeCoverage)
add_code_coverage_target(coverage make test)

# Add coverage flags (if debug build
if ( ENABLE_CODECOV )
   SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_COVERAGE}" )
   SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_COVERAGE}" )
endif()
