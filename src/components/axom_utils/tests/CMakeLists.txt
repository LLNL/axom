#
# Specify list of tests
#

set(gtest_utils_tests
   utils_TaskTimer.cpp
   utils_fileUtilities.cpp
   utils_types.cpp

   numerics_determinants.cpp
   numerics_linear_solve.cpp
   numerics_lu.cpp
   numerics_matrix.cpp
   )


set(utils_tests_depends
   axom_utils
   gtest
   )


#
# Add gtest based tests
#
foreach(test ${gtest_utils_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${utils_tests_depends} )
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name} )
endforeach()



#
# Add a test to query properties of the configuration
# Depends on boost to get its version
#

set(utils_config_test_depends ${utils_tests_depends})

if(BOOST_FOUND)
    list(APPEND utils_config_test_depends boost)
endif()

blt_add_executable( NAME utils_config
                    SOURCES utils_config.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${utils_config_test_depends}
                    )

blt_add_test( NAME utils_config
              COMMAND utils_config )

# Add flag to ignore unknown openmp pragmas in utils_config
set_target_properties(utils_config
                      PROPERTIES COMPILE_FLAGS ${AXOM_DISABLE_OMP_PRAGMA_WARNINGS})


if(ENABLE_OPENMP)
    # Add a config test that uses OpenMP
    blt_add_executable( NAME utils_config_openmp
                        SOURCES utils_config.cpp
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${utils_config_test_depends}
                        USE_OPENMP True)

    blt_add_test( NAME utils_config_openmp
                  COMMAND utils_config_openmp )
endif()

