#
# Specify list of tests
#

set(gtest_slic_tests
   testSlicAsserts.cpp
   )

set(bench_slic_tests
    benchSlicAsserts.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_slic_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON slic)
endforeach()


foreach(test ${bench_slic_tests})
  add_benchmark( TEST_SOURCE ${test} CLAS "--benchmark_min_time=0.00000001" DEPENDS_ON slic)
endforeach()
