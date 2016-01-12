#
# Specify list of tests
#

set(gtest_slic_tests
   slic_asserts.cpp
   )

set(bench_slic_tests
    slic_benchAsserts.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_slic_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON slic)
endforeach()


foreach(test ${bench_slic_tests})
  add_benchmark( TEST_SOURCE ${test} TEST_ARGS "--benchmark_min_time=0.00000001" DEPENDS_ON slic)
endforeach()
