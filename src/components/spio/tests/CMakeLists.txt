#
# Specify list of tests
#

make_executable(NAME spio_basicWriteRead
                SOURCES spio_basicWriteRead.cpp
                DEPENDS_ON spio 
                ADD_CTEST)
