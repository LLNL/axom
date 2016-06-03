#
# Specify list of tests
#

set(gtest_common_tests
   common_TaskTimer.cpp
   common_fileUtilities.cpp
   common_types.cpp
   )


#
# Add gtest based tests
#
foreach(test ${gtest_common_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON common gtest )
    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name} )
endforeach()



