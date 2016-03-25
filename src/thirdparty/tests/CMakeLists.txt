################################
# Add the conduit include path
################################
include_directories(${CONDUIT_INCLUDE_DIRS})

#
# Specify list of tests
#

################
# C++ Tests
################
set(gtest_thirdparty_tests
    gtest_smoke.cpp
    conduit_smoke.cpp
    )

set(gbenchmark_thirdparty_tests
    benchmark_smoke.cpp
    )

################
# Fortran Tests
################
set(gtest_thirdparty_fortran_tests
    fruit_smoke.f
    f_conduit_smoke.f)


################
# Add gtest based tests
################
foreach(test ${gtest_thirdparty_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON conduit)
endforeach()

if (ENABLE_BENCHMARK)
  foreach(test ${gbenchmark_thirdparty_tests})
      add_benchmark( TEST_SOURCE ${test} TEST_ARGS "--benchmark_min_time=0.0001")
  endforeach()
endif()


################
# Add fortran tests
################
if (ENABLE_FORTRAN)
  foreach(test ${gtest_thirdparty_fortran_tests})
      add_fortran_test(
          TEST_SOURCE ${test}
          DEPENDS_ON fruit conduit conduit_fortran
      )
  endforeach()
endif()


################
# Add OpenMP tests
################
if (ENABLE_OPENMP)
    make_executable( EXECUTABLE_NAME openmp_smoke 
                     EXECUTABLE_SOURCE openmp_smoke.cpp 
                     USE_OPENMP TRUE 
                     ADD_CTEST)
endif()


################
# Add MPI tests
################
if (ENABLE_MPI)
    add_mpi_gtest(TEST_SOURCE mpi_smoke.cpp NUM_PROCS 2)
endif()

