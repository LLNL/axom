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
    DEPENDS_ON sidre conduit conduit_relay hdf5 slic sparsehash
    )
blt_add_test( NAME sidre_lulesh2
              COMMAND sidre_lulesh2 -s 10
              )
