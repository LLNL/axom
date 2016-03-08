# 
# Build lulesh example
# NOTE: This should use a standard macro so we can enable features like C++11, MPI, OpenMP etc...
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

set_source_files_properties( ${lulesh_sources} PROPERTIES LANGUAGE "CXX" )

add_executable( sidre_lulesh2_example  ${lulesh_sources} ${lulesh_headers})
target_link_libraries( sidre_lulesh2_example  sidre  conduit slic )

add_test( NAME sidre_lulesh2_example
          COMMAND sidre_lulesh2_example -s 10
          WORKING_DIRECTORY ${TEST_OUTPUT_DIRECTORY}
        )
