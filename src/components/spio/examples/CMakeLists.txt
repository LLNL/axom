
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
    make_executable(
        NAME spio_${exe_name}
        SOURCES ${example_source}
        DEPENDS_ON spio sidre conduit conduit_io slic ${EXTRA_LIBS}
        IS_EXAMPLE)
endforeach()

