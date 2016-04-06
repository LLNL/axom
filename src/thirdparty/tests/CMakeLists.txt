################################
# Add the conduit include path
################################
include_directories(${CONDUIT_INCLUDE_DIRS})

#
# Specify list of tests
#

################
# Add gtest tests
################
add_gtest( TEST_SOURCE gtest_smoke.cpp )


################
# Test conduit tests
################
if (CONDUIT_FOUND)
  add_gtest( TEST_SOURCE conduit_smoke.cpp DEPENDS_ON conduit)

  if (ENABLE_FORTRAN)
    add_fortran_test(
        TEST_SOURCE f_conduit_smoke.f
        DEPENDS_ON fruit conduit conduit_fortran
    )
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
  make_executable( NAME openmp_smoke 
                   SOURCES openmp_smoke.cpp 
                   USE_OPENMP TRUE 
                   ADD_CTEST)
endif()


################
# Add MPI tests
################
if (ENABLE_MPI)
    add_mpi_gtest(TEST_SOURCE mpi_smoke.cpp NUM_PROCS 2)
endif()

