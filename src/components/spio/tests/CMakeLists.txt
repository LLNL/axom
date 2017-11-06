###############################################################################
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#
##############################################################################

#
# Specify list of tests
#

set(gtest_spio_serial_tests
   spio_serial.cpp
   )

set(gtest_spio_parallel_tests
   spio_parallel.cpp
   )

set(spio_test_depends
    spio
    conduit
    hdf5
    )

if(SCR_FOUND)
     list(APPEND spio_test_depends scr)
     list(APPEND gtest_spio_parallel_tests spio_scr.cpp)
endif()

if (ENABLE_MPI)

   foreach(test ${gtest_spio_serial_tests})
       get_filename_component( test_name ${test} NAME_WE )
       blt_add_executable( NAME ${test_name}_test
                           SOURCES ${test}
                           OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                           DEPENDS_ON ${spio_test_depends} gtest
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
                           )

       blt_add_test( NAME ${test_name} 
                     COMMAND ${test_name}_test
                     NUM_MPI_TASKS 4
                     )
   endforeach()

endif()

if(ENABLE_FORTRAN)
   if (ENABLE_MPI)
      blt_add_executable(NAME spio_basicWriteRead_F_test
                         SOURCES F_spio_basicWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         )
      blt_add_test(NAME F_spio_basicWriteRead
                   COMMAND spio_basicWriteRead_F_test
                   NUM_MPI_TASKS 1
                   )
                   
      blt_add_executable(NAME spio_parallelWriteRead_F_test
                         SOURCES F_spio_parallelWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         )
      blt_add_test(NAME F_spio_parallelWriteRead
                   COMMAND spio_parallelWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )
                   
      blt_add_executable(NAME spio_irregularWriteRead_F_test
                         SOURCES F_spio_irregularWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         )
      blt_add_test(NAME F_spio_irregularWriteRead
                   COMMAND spio_irregularWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )
                   
      blt_add_executable(NAME spio_externalWriteRead_F_test
                         SOURCES F_spio_externalWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         )
      blt_add_test(NAME F_spio_externalWriteRead
                   COMMAND spio_externalWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )

      blt_add_executable(NAME spio_preserveWriteRead_F_test
                         SOURCES F_spio_preserveWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends} fruit
                         )
      blt_add_test(NAME F_spio_preserveWriteRead
                   COMMAND spio_preserveWriteRead_F_test
                   NUM_MPI_TASKS 4
                   )

   endif()
endif()

