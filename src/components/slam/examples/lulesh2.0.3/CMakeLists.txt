######################################
# Slam (MeshAPI) version of lulesh 2.0
######################################


set(lulesh_headers
    lulesh.hpp
    # lulesh_tuple.hpp
    )

set(lulesh_sources
    lulesh.cpp
    lulesh-comm.cpp
    lulesh-init.cpp
    lulesh-util.cpp
    lulesh-viz.cpp
    )


# Ignore omp pragmas when openMP is not defined
if ( NOT ENABLE_OPENMP )
    set_source_files_properties(${lulesh_sources} 
        PROPERTIES COMPILE_FLAGS "${ATK_DISABLE_OMP_PRAGMA_WARNINGS}" )
endif()


blt_add_executable( 
    NAME        slam_lulesh
    SOURCES     ${lulesh_sources} ${lulesh_headers}
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  slic slam
    USE_OPENMP  ${ENABLE_OPENMP}
    )

blt_add_test(
    NAME        slam_lulesh
    COMMAND     slam_lulesh
    )
