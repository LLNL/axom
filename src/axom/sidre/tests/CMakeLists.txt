# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#
# Specify list of tests
#

set(gtest_sidre_tests
   sidre_smoke.cpp
   sidre_buffer.cpp
   sidre_buffer_unit.cpp
   sidre_class.cpp
   sidre_datastore.cpp
   sidre_datastore_unit.cpp
   sidre_external.cpp
   sidre_group.cpp
   sidre_opaque.cpp
   sidre_view.cpp
   sidre_native_layout.cpp
   sidre_attribute.cpp
   sidre_iterator.cpp
   sidre_array.cpp
   )

set(gtest_sidre_C_tests
   sidre_smoke_C.cpp
   sidre_buffer_C.cpp
   sidre_external_C.cpp
   sidre_group_C.cpp
   sidre_opaque_C.cpp
   sidre_types_C.cpp
   sidre_view_C.cpp
   )

set(fruit_sidre_tests
   sidre_smoke_F.f
   sidre_group_F.F
   sidre_buffer_F.f
   sidre_view_F.f
   sidre_external_F.f
   sidre_opaque_F.f
   sidre_allocatable_F.f
)


set(sidre_test_depends
    sidre
    conduit::conduit
    slic
    )

blt_list_append(TO sidre_test_depends ELEMENTS hdf5 IF ${HDF5_FOUND})

#
# Add gtest C++ based tests
#
foreach(test ${gtest_sidre_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}_test
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${sidre_test_depends} gtest
                        FOLDER axom/sidre/tests
                        )

    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}_test
                  )
endforeach()

#
# Add gtest C based tests
#
if (ENABLE_FORTRAN)
    foreach(test ${gtest_sidre_C_tests})
        get_filename_component( test_name ${test} NAME_WE )
        blt_add_executable( NAME ${test_name}_test
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON ${sidre_test_depends} gtest
                            FOLDER axom/sidre/tests
                            )
        blt_add_test( NAME ${test_name}
                      COMMAND ${test_name}_test
                      )
    endforeach()
endif()

#
# Add fortran tests
#
if (ENABLE_FORTRAN)
    foreach(test ${fruit_sidre_tests})
        get_filename_component( test_name ${test} NAME_WE )
        blt_add_executable( NAME ${test_name}_test
                            SOURCES ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON ${sidre_test_depends} fruit
                            FOLDER axom/sidre/tests
                            )
        blt_add_test( NAME ${test_name}
                      COMMAND ${test_name}_test
                      )
    endforeach()
endif()

if (ENABLE_MPI)
    add_subdirectory(spio)
endif()
