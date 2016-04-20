#
# Specify list of tests
#

blt_add_executable(NAME spio_basicWriteRead
                   SOURCES spio_basicWriteRead.cpp
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON spio conduit_io hdf5
                   )
blt_add_test(NAME spio_basicWriteRead
             COMMAND spio_basicWriteRead
             )
