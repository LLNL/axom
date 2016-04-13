#
# Specify list of tests
#

make_executable(ADD_CTEST
                EXECUTABLE_SOURCE spio_basicWriteRead.cpp
                DEPENDS_ON spio )
if (ENABLE_MPI)
    add_mpi_gtest(TEST_SOURCE spio_parallelWriteRead.cpp
                  NUM_PROCS 4
                  DEPENDS_ON spio )
endif()

