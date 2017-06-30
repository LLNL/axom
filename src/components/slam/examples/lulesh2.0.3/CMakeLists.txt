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
        PROPERTIES COMPILE_FLAGS "${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}" )
endif()


# blt_add_executable depends_on handles our include dir settings
# these teste are using boost so boost needs to be a dep
set(slam_tests_depends_on axom_utils slic slam)
if (BOOST_FOUND)
    list(APPEND slam_tests_depends_on boost)
endif()

if (ENABLE_OPENMP)
    list(APPEND slam_tests_depends_on openmp)
endif()

blt_add_executable(
    NAME        slam_lulesh
    SOURCES     ${lulesh_sources} ${lulesh_headers}
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${slam_tests_depends_on}
    )

if(ENABLE_MPI)
    blt_add_test(NAME        slam_lulesh
                 COMMAND     slam_lulesh
                 NUM_PROCS   8
                )
else()
    blt_add_test(NAME        slam_lulesh
                 COMMAND     slam_lulesh
                )
endif()

