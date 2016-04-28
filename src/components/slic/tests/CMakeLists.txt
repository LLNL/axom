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
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON slic gtest )
    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name} )
endforeach()


foreach(test ${bench_slic_tests})
	blt_add_benchmark( TEST_SOURCE ${test} TEST_ARGS "--benchmark_min_time=0.00000001" DEPENDS_ON slic )
endforeach()
