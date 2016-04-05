####################################
# Original version of lulesh 2.0
# with a few light edits
####################################


set(lulesh_lib_headers
    lulesh.h
    #lulesh_tuple.hpp
    )

set(lulesh_lib_sources
    lulesh-comm.cc
    lulesh-init.cc
    lulesh-util.cc
    lulesh-viz.cc
    )
set(lulesh_exe_source
    lulesh.cc
    )

# Add some extra configuration flags for MPI, OpenMP and allowing unused parameters (in 3rd party source code)
set(luleshExtraCompileFlags ${ATK_DISABLE_UNUSED_PARAMETER_WARNINGS})
if ( NOT ENABLE_OPENMP )
    set(luleshExtraCompileFlags "${luleshExtraCompileFlags} ${ATK_DISABLE_OMP_PRAGMA_WARNINGS}")
endif()

if ( ENABLE_MPI )
    set(luleshExtraCompileFlags "${luleshExtraCompileFlags} -DUSE_MPI=1 ")
else()
    set(luleshExtraCompileFlags "${luleshExtraCompileFlags} -DUSE_MPI=0 ")    
endif()
#message(STATUS " luleshExtraCompileFlags is '${luleshExtraCompileFlags}'")

set_source_files_properties(${lulesh_lib_sources} ${lulesh_lib_headers} ${lulesh_exe_source}  
                            PROPERTIES COMPILE_FLAGS " ${luleshExtraCompileFlags} " )

make_library( 
    LIBRARY_NAME slamLuleshOrigLib 
    LIBRARY_SOURCES "${lulesh_lib_sources}"  "${lulesh_lib_headers}"  
    USE_OPENMP ${ENABLE_OPENMP} 
   ) 
   
make_executable( 
    EXECUTABLE_NAME slam_OrigLulesh 
    EXECUTABLE_SOURCE ${lulesh_exe_source} 
    DEPENDS_ON "slamLuleshOrigLib" 
    USE_OPENMP ${ENABLE_OPENMP}
    IS_EXAMPLE
   ) 
      
