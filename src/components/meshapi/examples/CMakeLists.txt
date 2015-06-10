
#
# List of example sources
#
set(example_sources
	meshapiShockTube.cpp
    UnstructMeshField.cpp
    )



#
# Build the example sources
#
foreach(example_source ${example_sources})
    make_executable(EXECUTABLE_SOURCE ${example_source} DEPENDS_ON meshapi )
endforeach()


# 
# Build lulesh example
# NOTE: This should use a standard macro so we can enable features like C++11, MPI, OpenMP etc...
#

set(lulesh_dir lulesh2.0.3  )

set(lulesh_headers
    ${lulesh_dir}/lulesh.hpp
    # ${lulesh_dir}/lulesh_tuple.hpp
    )

set(lulesh_sources
    ${lulesh_dir}/lulesh.cpp
    ${lulesh_dir}/lulesh-comm.cpp
    ${lulesh_dir}/lulesh-init.cpp
    ${lulesh_dir}/lulesh-util.cpp
    ${lulesh_dir}/lulesh-viz.cpp
    )

## Ignore omp pragmas when openMP is not defined
set_source_files_properties(${lulesh_sources} PROPERTIES COMPILE_FLAGS " -Wno-unknown-pragmas " )
add_executable( meshapiLulesh  ${lulesh_sources} ${lulesh_headers})
target_link_libraries( meshapiLulesh  meshapi )