
#
# List of example sources
#
set(example_sources
    ShockTube.cpp
    DumpStateFile.cpp
    )

#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable( EXECUTABLE_SOURCE ${example_source} DEPENDS_ON sidre conduit slic )
endforeach()

# 
# Build lulesh example
# NOTE: This should use a standard macro so we can enable features like C++11, MPI, OpenMP etc...
#

set(lulesh_dir lulesh2  )

set(lulesh_headers
    ${lulesh_dir}/lulesh.h
    ${lulesh_dir}/lulesh_tuple.h
    )

set(lulesh_sources
    ${lulesh_dir}/lulesh.cc
    ${lulesh_dir}/lulesh-comm.cc
    ${lulesh_dir}/lulesh-init.cc
    ${lulesh_dir}/lulesh-util.cc
    ${lulesh_dir}/lulesh-viz.cc
    )

set_source_files_properties( ${lulesh_sources} PROPERTIES LANGUAGE "CXX" )

add_executable( sidreLulesh2  ${lulesh_sources} ${lulesh_headers})
target_link_libraries( sidreLulesh2  sidre  conduit slic )
