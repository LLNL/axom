####################################
# Original version of lulesh 2.0
# with a few light edits
####################################


set(lulesh_headers
        lulesh.h
        #lulesh_tuple.hpp
        )

set(lulesh_sources
        lulesh.cc
        lulesh-comm.cc
        lulesh-init.cc
        lulesh-util.cc
        lulesh-viz.cc
        )

# Add some extra configuration flags for MPI, OpenMP and allowing unused parameters (in 3rd party source code)
set(luleshExtraCompileFlags ${AXOM_DISABLE_UNUSED_PARAMETER_WARNINGS})
if ( NOT ENABLE_OPENMP )
    set(luleshExtraCompileFlags "${luleshExtraCompileFlags} ${AXOM_DISABLE_OMP_PRAGMA_WARNINGS}")
endif()

if ( ENABLE_MPI )
    set(luleshExtraCompileFlags "${luleshExtraCompileFlags} -DUSE_MPI=1 ")
else()
    set(luleshExtraCompileFlags "${luleshExtraCompileFlags} -DUSE_MPI=0 ")
endif()
#message(STATUS " luleshExtraCompileFlags is '${luleshExtraCompileFlags}'")

set_source_files_properties(${lulesh_sources} ${lulesh_headers} 
    PROPERTIES COMPILE_FLAGS " ${luleshExtraCompileFlags} " )


blt_add_executable(
    NAME        slam_lulesh_orig_ex
    SOURCES     ${lulesh_sources} ${lulesh_headers}
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    USE_OPENMP  ${ENABLE_OPENMP}
    ) 
