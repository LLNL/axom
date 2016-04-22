####################################
# Slam (MeshAPI) version of lulesh 2.0
####################################


set(lulesh_lib_headers
    lulesh.hpp
    # lulesh_tuple.hpp
    )

set(lulesh_lib_sources
    lulesh-comm.cpp
    lulesh-init.cpp
    lulesh-util.cpp
    lulesh-viz.cpp
    )

set(lulesh_exe_source
    lulesh.cpp
    )


## Ignore omp pragmas when openMP is not defined
if ( NOT ENABLE_OPENMP )
    set_source_files_properties(${lulesh_lib_sources} ${lulesh_exe_source}
        PROPERTIES COMPILE_FLAGS "${ATK_DISABLE_OMP_PRAGMA_WARNINGS}" )
endif()

blt_add_library( NAME      slamLuleshLib
                 SOURCES   "${lulesh_lib_sources}" "${lulesh_lib_headers}"
                 DEPENDS_ON slam
                 USE_OPENMP ${ENABLE_OPENMP}
                 ) 

blt_add_executable( NAME       slam_Lulesh
                     SOURCES    ${lulesh_exe_source}
                     OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                     DEPENDS_ON
                        slic
                        slam
                        slamLuleshLib
                      USE_OPENMP ${ENABLE_OPENMP}
                      )
blt_add_test(NAME slam_Lulesh
             COMMAND slam_Lulesh
             )
