#################################################
# Setup code metrics - coverage, profiling, etc
#################################################

################################
# Enable code coverage via gcov
# Note: Only supported for gnu.
################################

# switch to this one
#if(ENABLE_CODECOV)
#   include(CodeCoverage)
#endif()

# These should be set in some separate macro later that sets coverage flags for each compiler.
#SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
#SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
SET(GCC_COVERAGE_COMPILE_FLAGS "--coverage")
SET(GCC_COVERAGE_LINK_FLAGS    "--coverage")

if (CMAKE_BUILD_TYPE MATCHES Debug)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        MESSAGE(STATUS "Debug build is using gnu, code coverage via gcov enabled.")
        SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
        SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )
    else()
        MESSAGE(STATUS "Debug build is not using gnu, code coverage is disabled.")
        SET(ENABLE_CODECOV OFF)
    endif()
endif()
