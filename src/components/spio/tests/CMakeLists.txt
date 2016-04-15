#
# Specify list of tests
#

make_executable(ADD_CTEST
                EXECUTABLE_SOURCE spio_basicWriteRead.cpp
                DEPENDS_ON spio )
