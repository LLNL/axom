#
# Specify list of tests
#

set(spio_test_depends
    spio
    conduit
    hdf5
    )

if(SCR_FOUND)
     set(spio_test_depends
         ${spio_test_depends}
         scr
         )
endif()

if (ENABLE_MPI)
   blt_add_executable(NAME spio_basicWriteRead_test
                      SOURCES spio_basicWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${spio_test_depends}
                      )
   blt_add_test(NAME spio_basicWriteRead
                COMMAND spio_basicWriteRead_test
                NUM_PROCS 1
                )

   blt_add_executable(NAME spio_parallelWriteRead_test
                      SOURCES spio_parallelWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${spio_test_depends}
                      )
   blt_add_test(NAME spio_parallelWriteRead
                COMMAND spio_parallelWriteRead_test
                NUM_PROCS 4
                )

   blt_add_executable(NAME spio_SCRWriteRead_test
                      SOURCES spio_SCRWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${spio_test_depends}
                      )
   blt_add_test(NAME spio_SCRWriteRead
                COMMAND spio_SCRWriteRead_test
                NUM_PROCS 4
                )
                
   blt_add_executable(NAME spio_irregularWriteRead_test
                      SOURCES spio_irregularWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${spio_test_depends}
                      )
   blt_add_test(NAME spio_irregularWriteRead
                COMMAND spio_irregularWriteRead_test
                NUM_PROCS 4
                )
                
   blt_add_executable(NAME spio_externalWriteRead_test
                      SOURCES spio_externalWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${spio_test_depends}
                      )
   blt_add_test(NAME spio_externalWriteRead
                COMMAND spio_externalWriteRead_test
                NUM_PROCS 4
                )
                
   blt_add_executable(NAME spio_WRW_test
                      SOURCES spio_parallel_WRW.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON spio conduit hdf5
                      )
   blt_add_test(NAME spio_WRW_1
                COMMAND spio_WRW_test
                NUM_PROCS 1
                )
     
   blt_add_test(NAME spio_WRW_4
                COMMAND spio_WRW_test
                NUM_PROCS 4
                )

   blt_add_executable(NAME spio_preserveWriteRead_test
                      SOURCES spio_preserveWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${spio_test_depends}
                      )
   blt_add_test(NAME spio_preserveWriteRead
                COMMAND spio_preserveWriteRead_test
                NUM_PROCS 4
                )

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
                   NUM_PROCS 1
                   )
                   
      blt_add_executable(NAME spio_parallelWriteRead_F_test
                         SOURCES F_spio_parallelWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         )
      blt_add_test(NAME F_spio_parallelWriteRead
                   COMMAND spio_parallelWriteRead_F_test
                   NUM_PROCS 4
                   )
                   
      blt_add_executable(NAME spio_irregularWriteRead_F_test
                         SOURCES F_spio_irregularWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         )
      blt_add_test(NAME F_spio_irregularWriteRead
                   COMMAND spio_irregularWriteRead_F_test
                   NUM_PROCS 4
                   )
                   
      blt_add_executable(NAME spio_externalWriteRead_F_test
                         SOURCES F_spio_externalWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends}
                         )
      blt_add_test(NAME F_spio_externalWriteRead
                   COMMAND spio_externalWriteRead_F_test
                   NUM_PROCS 4
                   )

      blt_add_executable(NAME spio_preserveWriteRead_F_test
                         SOURCES F_spio_preserveWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON ${spio_test_depends} fruit
                         )
      blt_add_test(NAME F_spio_preserveWriteRead
                   COMMAND spio_preserveWriteRead_F_test
                   NUM_PROCS 4
                   )

   endif()
endif()

