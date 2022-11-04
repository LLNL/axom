# Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Spio unit tests
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_spio_serial_tests
   spio_serial.hpp
   )

set(gtest_spio_parallel_tests
   spio_basic.hpp
   spio_parallel.hpp
   )

set(spio_tests_depends_on axom gtest)
blt_list_append( TO spio_tests_depends_on ELEMENTS cuda IF ${ENABLE_CUDA} )

#------------------------------------------------------------------------------
# Add gtest C++ tests
#------------------------------------------------------------------------------
blt_add_executable( NAME       spio_tests
                    SOURCES    spio_main.cpp
                    HEADERS    ${gtest_spio_parallel_tests} ${gtest_spio_serial_tests}
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${spio_tests_depends_on}
                    FOLDER     axom/sidre/tests
                    )

foreach(test ${gtest_spio_serial_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_test( NAME          ${test_name} 
                   COMMAND       spio_tests --gtest_filter=${test_name}*
                   NUM_MPI_TASKS 1
                   )
endforeach()

foreach(test ${gtest_spio_parallel_tests})
    get_filename_component( test_name ${test} NAME_WE )
    axom_add_test( NAME          ${test_name} 
                   COMMAND       spio_tests --gtest_filter=${test_name}*
                   NUM_MPI_TASKS 4
                   )
endforeach()

#------------------------------------------------------------------------------
# Add Fortran tests
#------------------------------------------------------------------------------
set(spio_fortran_tests_depends_on sidre)
blt_list_append( TO spio_fortran_tests_depends_on ELEMENTS cuda IF ${ENABLE_CUDA} )

if(ENABLE_FORTRAN)
    # Non-Fruit parallel tests with 1 MPI task
    set(spio_one_fortran_tests
        F_spio_basicWriteRead.F
        F_spio_blueprintIndex.F
        )

    foreach(test ${spio_one_fortran_tests})
        get_filename_component( test_name ${test} NAME_WE )
        string(SUBSTRING "${test_name}" 2 -1 test_name) # remove F_
        blt_add_executable( NAME       ${test_name}_F_test
                            SOURCES    ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON ${spio_fortran_tests_depends_on}
                            FOLDER     axom/sidre/tests
                            )
        axom_add_test( NAME    ${test_name}_F
                      COMMAND ${test_name}_F_test
                      NUM_MPI_TASKS 1
                      )
    endforeach()

    # Non-Fruit parallel tests with 4 MPI task
    set(spio_four_fortran_tests
        F_spio_parallelWriteRead.F
        F_spio_irregularWriteRead.F
        F_spio_externalWriteRead.F
        )

    foreach(test ${spio_four_fortran_tests})
        get_filename_component( test_name ${test} NAME_WE )
        string(SUBSTRING "${test_name}" 2 -1 test_name) # remove F_
        blt_add_executable( NAME       ${test_name}_F_test
                            SOURCES    ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON ${spio_fortran_tests_depends_on}
                            FOLDER     axom/sidre/tests
                            )
        axom_add_test( NAME    ${test_name}_F
                       COMMAND ${test_name}_F_test
                       NUM_MPI_TASKS 4
                       )
    endforeach()

    # Fruit parallel tests with 4 MPI task
    set(spio_fruit_tests
        F_spio_preserveWriteRead.F
        )

    foreach(test ${spio_fruit_tests})
        get_filename_component( test_name ${test} NAME_WE )
        string(SUBSTRING "${test_name}" 2 -1 test_name) # remove F_
        blt_add_executable( NAME       ${test_name}_F_test
                            SOURCES    ${test}
                            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                            DEPENDS_ON ${spio_fortran_tests_depends_on} fruit
                            FOLDER     axom/sidre/tests
                            )
        axom_add_test( NAME    ${test_name}_F
                       COMMAND ${test_name}_F_test
                       NUM_MPI_TASKS 4
                       )
    endforeach()

endif()
