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
        add_fortran_test(
            TEST_SOURCE f_conduit_smoke.f
            DEPENDS_ON fruit conduit conduit_fortran)
  endif()
endif()


################
# Test gbenchmark
################
if (ENABLE_BENCHMARKS)
    add_benchmark( TEST_SOURCE benchmark_smoke.cpp TEST_ARGS "--benchmark_min_time=0.0001")
endif()


################
# Add fruit tests
################
if (ENABLE_FORTRAN)
    add_fortran_test(
      TEST_SOURCE fruit_smoke.f
      DEPENDS_ON fruit
  )
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

