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

if(ENABLE_BENCHMARK)
  foreach(test ${gbenchmark_thirdparty_tests})
      add_benchmark( TEST_SOURCE ${test} CLAS "--benchmark_min_time=0.01" )
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
# This is a special case, I am not sure we have logic to add WITH_OPENMP to gtest?
################
if(ENABLE_OMP)
    make_executable( EXECUTABLE_NAME openmp_smoke EXECUTABLE_SOURCE openmp_smoke.cpp WITH_OPENMP)
endif()



    

