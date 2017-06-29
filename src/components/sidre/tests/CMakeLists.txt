#
# Specify list of tests
#

set(gtest_sidre_tests
   sidre_smoke_test.cpp
   sidre_buffer_test.cpp
   sidre_buffer_unit_test.cpp
   sidre_class_test.cpp
   sidre_datastore_test.cpp
   sidre_datastore_unit_test.cpp
   sidre_external_test.cpp
   sidre_group_test.cpp
   sidre_opaque_test.cpp
   sidre_view_test.cpp
   sidre_native_layout_test.cpp
   sidre_attribute_test.cpp
   )

set(gtest_sidre_C_tests
   sidre_smoke_C_test.cpp
   sidre_buffer_C_test.cpp
   sidre_external_C_test.cpp
   sidre_group_C_test.cpp
   sidre_opaque_C_test.cpp
   sidre_types_C_test.cpp
   sidre_view_C_test.cpp
   )

set(fruit_sidre_tests
   sidre_smoke_F_test.f
   sidre_group_F_test.f
   sidre_buffer_F_test.f
   sidre_view_F_test.f
   sidre_external_F_test.f
   sidre_opaque_F_test.f
   sidre_allocatable_F_test.f
)

#
# Add gtest C++ based tests
#
foreach(test ${gtest_sidre_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON sidre conduit conduit_relay hdf5 slic gtest
                        )

    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name}
                  )
endforeach()

#
# Add gtest C based tests
#
if (ENABLE_FORTRAN)
    foreach(test ${gtest_sidre_C_tests})
        get_filename_component( test_name ${test} NAME_WE )
        blt_add_executable( NAME ${test_name}
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON sidre conduit conduit_relay hdf5 slic gtest
                            )
        blt_add_test( NAME ${test_name} 
                      COMMAND ${test_name}
                      )
    endforeach()
endif()

#
# Add fortran tests
#
if (ENABLE_FORTRAN)
    foreach(test ${fruit_sidre_tests})
        get_filename_component( test_name ${test} NAME_WE )
        blt_add_executable( NAME ${test_name}
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON sidre conduit conduit_relay hdf5 slic fruit
                            )
        blt_add_test( NAME ${test_name} 
                      COMMAND ${test_name}
                      )
    endforeach()
endif()

