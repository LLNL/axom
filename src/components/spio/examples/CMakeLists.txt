
#
# List of single source file examples
#
set(example_sources
    IORead.cpp
    IOWrite.cpp
)



#
# Build the example sources
#
foreach(example_source ${example_sources})
    get_filename_component(exe_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME spio_${exe_name}
        SOURCES ${example_source}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON spio sidre conduit conduit_relay hdf5 slic ${EXTRA_LIBS}
        )
endforeach()

