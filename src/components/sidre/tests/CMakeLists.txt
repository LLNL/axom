#
# Specify list of tests
#

set(gtest_sidre_tests
   sidre_smoke.cpp
   sidre_buffer.cpp
   sidre_class.cpp
   sidre_datastore.cpp
   sidre_external.cpp
   sidre_group.cpp
   sidre_opaque.cpp
   sidre_view.cpp
   sidre_native_layout.cpp
   )

set(gtest_sidre_C_tests
   C_sidre_smoke.cpp
   C_sidre_buffer.cpp
   C_sidre_external.cpp
   C_sidre_group.cpp
   C_sidre_opaque.cpp
   C_sidre_types.cpp
   C_sidre_view.cpp
   )

set(fruit_sidre_tests
   F_sidre_smoke.f
   F_sidre_group.f
   F_sidre_buffer.f
   F_sidre_view.f
   F_sidre_external.f
   F_sidre_opaque.f
   F_sidre_allocatable.f
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

