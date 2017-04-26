#
# Specify list of tests
#

set(gtest_common_tests
   common_TaskTimer.cpp
   common_fileUtilities.cpp
   common_types.cpp

   numerics_determinants.cpp
   numerics_lu.cpp
   numerics_matrix.cpp
   )


set(common_tests_depends
   axom_utils
   gtest
   )


#
# Add gtest based tests
#
foreach(test ${gtest_common_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${common_tests_depends} )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name} )
endforeach()



#
# Add a test to query properties of the configuration
# Depends on boost to get its version
#

set(common_config_test_depends ${common_tests_depends})

if(BOOST_FOUND)
    list(APPEND common_config_test_depends boost)
endif()

blt_add_executable( NAME common_config
                    SOURCES common_config.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${common_config_test_depends}
                    )

blt_add_test( NAME common_config
              COMMAND common_config )

# Add flag to ignore unknown openmp pragmas in common_config
set_target_properties(common_config
                      PROPERTIES COMPILE_FLAGS ${AXOM_DISABLE_OMP_PRAGMA_WARNINGS})


if(ENABLE_OPENMP)
    # Add a config test that uses OpenMP
    blt_add_executable( NAME common_config_openmp
                        SOURCES common_config.cpp
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${common_config_test_depends}
                        USE_OPENMP True)

    blt_add_test( NAME common_config_openmp
                  COMMAND common_config_openmp )
endif()

