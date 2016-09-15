#
# Specify list of tests
#

blt_add_executable(NAME spio_basicWriteRead
                   SOURCES spio_basicWriteRead.cpp
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON spio conduit hdf5
                   )
blt_add_test(NAME spio_basicWriteRead
             COMMAND spio_basicWriteRead
             )

if (ENABLE_MPI)
   blt_add_executable(NAME spio_parallelWriteRead
                      SOURCES spio_parallelWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON spio conduit hdf5
                      )
   blt_add_test(NAME spio_parallelWriteRead
                COMMAND spio_parallelWriteRead
                NUM_PROCS 4
                )
   blt_add_executable(NAME spio_irregularWriteRead
                      SOURCES spio_irregularWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON spio conduit hdf5
                      )
   blt_add_test(NAME spio_irregularWriteRead
                COMMAND spio_irregularWriteRead
                NUM_PROCS 4
                )
   blt_add_executable(NAME spio_externalWriteRead
                      SOURCES spio_externalWriteRead.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON spio conduit hdf5
                      )
   blt_add_test(NAME spio_externalWriteRead
                COMMAND spio_externalWriteRead
                NUM_PROCS 4
                )
endif()

if(ENABLE_FORTRAN)
   blt_add_executable(NAME F_spio_basicWriteRead
                      SOURCES F_spio_basicWriteRead.f
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON spio conduit hdf5
                      )
   blt_add_test(NAME F_spio_basicWriteRead
                COMMAND F_spio_basicWriteRead
                )

   if (ENABLE_MPI)
      blt_add_executable(NAME F_spio_parallelWriteRead
                         SOURCES F_spio_parallelWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON spio conduit hdf5
                         )
      blt_add_test(NAME F_spio_parallelWriteRead
                   COMMAND F_spio_parallelWriteRead
                   NUM_PROCS 4
                   )
      blt_add_executable(NAME F_spio_irregularWriteRead
                         SOURCES F_spio_irregularWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON spio conduit hdf5
                         )
      blt_add_test(NAME F_spio_irregularWriteRead
                   COMMAND F_spio_irregularWriteRead
                   NUM_PROCS 4
                   )
      blt_add_executable(NAME F_spio_externalWriteRead
                         SOURCES F_spio_externalWriteRead.f
                         OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                         DEPENDS_ON spio conduit hdf5
                         )
      blt_add_test(NAME F_spio_externalWriteRead
                   COMMAND F_spio_externalWriteRead
                   NUM_PROCS 4
                   )
   endif()
endif()
