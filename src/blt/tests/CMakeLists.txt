#
# Smoke tests for BLT tpls
#

################
# gtest test
################
blt_add_executable(NAME blt_gtest_smoke
                   SOURCES blt_gtest_smoke.cpp
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON gtest)

blt_add_test( NAME blt_gtest_smoke
              COMMAND blt_gtest_smoke)


################
# gbenchmark test
################
if (ENABLE_BENCHMARKS)
    blt_add_benchmark( TEST_SOURCE blt_benchmark_smoke.cpp TEST_ARGS "--benchmark_min_time=0.0001")
endif()


################
# fruit test
################
if (ENABLE_FORTRAN)
    blt_add_executable(NAME blt_fruit_smoke
                       SOURCES blt_fruit_smoke.f
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       DEPENDS_ON fruit)

    blt_add_test(NAME blt_fruit_smoke
                 COMMAND blt_fruit_smoke)
endif()


################
# OpenMP test
################
if (ENABLE_OPENMP)
    blt_add_executable(NAME blt_openmp_smoke 
                       SOURCES blt_openmp_smoke.cpp 
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                       USE_OPENMP TRUE )

    blt_add_test(NAME blt_openmp_smoke
                 COMMAND blt_openmp_smoke
                 NUM_PROCS 2)
endif()


################
#  MPI test
################
if (ENABLE_MPI)
    blt_add_executable(NAME blt_mpi_smoke
                       SOURCES blt_mpi_smoke.cpp
                       OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY} )

    blt_add_test(NAME blt_mpi_smoke
                 COMMAND blt_mpi_smoke
                 NUM_PROCS 2)
endif()

