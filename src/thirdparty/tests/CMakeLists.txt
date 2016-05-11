#
# Specify list of tests
#

################
# Add gtest tests
################
blt_add_executable(NAME gtest_smoke
                   SOURCES gtest_smoke.cpp
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON gtest)
blt_add_test( NAME gtest_smoke
              COMMAND gtest_smoke)


################
# Test conduit tests
################
if (CONDUIT_FOUND)
    blt_add_executable(NAME conduit_smoke
                       SOURCES conduit_smoke.cpp
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON conduit gtest)
    blt_add_test(NAME conduit_smoke
                 COMMAND conduit_smoke)

    if (ENABLE_FORTRAN)
        blt_add_executable(NAME f_conduit_smoke
                           SOURCES f_conduit_smoke.f
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON conduit conduit_fortran fruit)
        blt_add_test(NAME f_conduit_smoke
                     COMMAND f_conduit_smoke)
  endif()
endif()


################
# Test gbenchmark
################
if (ENABLE_BENCHMARKS)
    blt_add_benchmark( TEST_SOURCE benchmark_smoke.cpp TEST_ARGS "--benchmark_min_time=0.0001")
endif()


################
# Add fruit tests
################
if (ENABLE_FORTRAN)
    blt_add_executable(NAME fruit_smoke
                       SOURCES fruit_smoke.f
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON fruit)
    blt_add_test(NAME fruit_smoke
                 COMMAND fruit_smoke)
endif()


################
# Add OpenMP tests
################
if (ENABLE_OPENMP)
    blt_add_executable(NAME openmp_smoke 
                       SOURCES openmp_smoke.cpp 
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       USE_OPENMP TRUE )
    blt_add_test(NAME openmp_smoke
                 COMMAND openmp_smoke
                 NUM_PROCS 2)
endif()


################
# Add MPI tests
################
if (ENABLE_MPI)
    blt_add_executable(NAME mpi_smoke
                       SOURCES mpi_smoke.cpp
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY} )
    blt_add_test(NAME mpi_smoke
                 COMMAND mpi_smoke
                 NUM_PROCS 2)
endif()

################################################################################
# BLT compiler flag tests
# Tests should fail if we have the wrong warning disabler set for the
# given compiler/configuration.  Flags are defined in blt's SetupCompilerOptions.
################################################################################

# Test for ATK_DISABLE_UNINITIALIZED_WARNINGS flag
set( BLT_SMOKE_TESTS
        blt_uninitialized.cpp
        blt_unused_param.cpp
        blt_unused_var.cpp
        blt_nopragma.cpp
        blt_strict_aliasing.cpp
        )

foreach(test ${BLT_SMOKE_TESTS})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( 
        NAME ${test_name}
        SOURCES ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        USE_OPENMP False
        )
    blt_add_test( 
        NAME ${test_name} 
        COMMAND ${test_name} 
        )
endforeach()

set( BLT_SMOKE_FLAGS "${ATK_TREAT_WARNINGS_AS_ERRORS} ${ATK_ENABLE_ALL_WARNINGS}")

set_target_properties(blt_uninitialized 
    PROPERTIES COMPILE_FLAGS "${BLT_SMOKE_FLAGS} ${ATK_DISABLE_UNINITIALIZED_WARNINGS}")
set_target_properties(blt_unused_param 
    PROPERTIES COMPILE_FLAGS "${BLT_SMOKE_FLAGS} ${ATK_DISABLE_UNUSED_PARAMETER_WARNINGS}")
set_target_properties(blt_unused_var 
    PROPERTIES COMPILE_FLAGS "${BLT_SMOKE_FLAGS} ${ATK_DISABLE_UNUSED_VARIABLE_WARNINGS}")
set_target_properties(blt_nopragma 
    PROPERTIES COMPILE_FLAGS "${BLT_SMOKE_FLAGS} ${ATK_DISABLE_OMP_PRAGMA_WARNINGS}")
    
# the aliasing test requires optimization, set this up as custom compiler flag
# in case we need more flexibility for this flag (e.g. for MSVC)
blt_append_custom_compiler_flag(FLAGS_VAR TMP_OPT_FLAG DEFAULT "-O2")
set_target_properties(blt_strict_aliasing 
    PROPERTIES COMPILE_FLAGS "${BLT_SMOKE_FLAGS} ${TMP_OPT_FLAG} ${ATK_DISABLE_ALIASING_WARNINGS}")


