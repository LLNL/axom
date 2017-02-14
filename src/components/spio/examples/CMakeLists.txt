
#
# List of single source file examples
#
set(example_sources
    IORead.cpp
    IOWrite.cpp
)

set(SPIO_EXAMPLE_DEPENDS spio sidre conduit conduit_relay hdf5 slic ${EXTRA_LIBS})

#
# Build the example sources
#
foreach(src ${example_sources})
    get_filename_component(exe_name ${src} NAME_WE)
    blt_add_executable(
        NAME spio_${exe_name}
        SOURCES ${src}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON ${SPIO_EXAMPLE_DEPENDS}
        )
endforeach()

