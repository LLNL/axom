#
# Build lulesh example
#

set(lulesh_headers
    lulesh.h
    lulesh_tuple.h
    )

set(lulesh_sources
    lulesh.cc
    lulesh-comm.cc
    lulesh-init.cc
    lulesh-util.cc
    lulesh-viz.cc
    )

blt_add_executable(
    NAME sidre_lulesh2
    SOURCES ${lulesh_sources} ${lulesh_headers}
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON axom_common sidre conduit conduit_relay hdf5 slic
    )

if(ENABLE_MPI)
    blt_add_test(NAME sidre_lulesh2
                 COMMAND sidre_lulesh2 -s 5
                 NUM_PROCS 8
                )
else()
    blt_add_test(NAME sidre_lulesh2
                 COMMAND sidre_lulesh2 -s 5
                )
endif()
