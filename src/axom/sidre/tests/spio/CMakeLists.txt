# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Spio unit tests
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_spio_serial_tests
   spio_serial.cpp
   )

set(gtest_spio_parallel_tests
   spio_basic.cpp
   spio_parallel.cpp
   )

set(spio_test_depends
    sidre
    conduit::conduit
    mpi
    )

blt_list_append(TO spio_test_depends ELEMENTS hdf5 IF HDF5_FOUND)

if(SCR_FOUND)
     list(APPEND spio_test_depends scr)
     list(APPEND gtest_spio_parallel_tests spio_scr.cpp)
endif()

#------------------------------------------------------------------------------
# Add gtest C++ tests
#------------------------------------------------------------------------------
if (ENABLE_MPI)

   foreach(test ${gtest_spio_serial_tests})
       get_filename_component( test_name ${test} NAME_WE )
       blt_add_executable( NAME ${test_name}_test
                           SOURCES ${test}
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON ${spio_test_depends} gtest
                           FOLDER axom/sidre/tests
                           )

       blt_add_test( NAME ${test_name} 
                     COMMAND ${test_name}_test
                     NUM_MPI_TASKS 1
                    )
   endforeach()


   foreach(test ${gtest_spio_parallel_tests})
       get_filename_component( test_name ${test} NAME_WE )
       blt_add_executable( NAME ${test_name}_test
                           SOURCES ${test}
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON ${spio_test_depends} gtest
                           FOLDER axom/sidre/tests
                           )

       blt_add_test( NAME ${test_name} 
                     COMMAND ${test_name}_test
                     NUM_MPI_TASKS 4
                     )
   endforeach()

endif()

#------------------------------------------------------------------------------
# Add Fortran tests
#------------------------------------------------------------------------------
if(ENABLE_FORTRAN)
   if (ENABLE_MPI)
      blt_add_executable(NAME spio_basicWriteRead_F_test
                         SOURCES F_spio_basicWriteRead.F
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         FOLDER axom/sidre/tests
                         )
      blt_add_test(NAME F_spio_basicWriteRead
                   COMMAND spio_basicWriteRead_F_test
                   NUM_MPI_TASKS 1
                   )

      blt_add_executable(NAME spio_parallelWriteRead_F_test
                         SOURCES F_spio_parallelWriteRead.F
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         FOLDER axom/sidre/tests
                         )
      blt_add_test(NAME F_spio_parallelWriteRead
                   COMMAND spio_parallelWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )

      blt_add_executable(NAME spio_irregularWriteRead_F_test
                         SOURCES F_spio_irregularWriteRead.F
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         FOLDER axom/sidre/tests
                         )
      blt_add_test(NAME F_spio_irregularWriteRead
                   COMMAND spio_irregularWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )

      blt_add_executable(NAME spio_externalWriteRead_F_test
                         SOURCES F_spio_externalWriteRead.F
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         FOLDER axom/sidre/tests
                         )
      blt_add_test(NAME F_spio_externalWriteRead
                   COMMAND spio_externalWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )

      blt_add_executable(NAME spio_preserveWriteRead_F_test
                         SOURCES F_spio_preserveWriteRead.F
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends} fruit
                         FOLDER axom/sidre/tests
                         )
      blt_add_test(NAME F_spio_preserveWriteRead
                   COMMAND spio_preserveWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )

   endif()
endif()

