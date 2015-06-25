
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

set(lulesh_lib_headers
    ${lulesh_dir}/lulesh.hpp
    # ${lulesh_dir}/lulesh_tuple.hpp
    )

set(lulesh_lib_sources
    ${lulesh_dir}/lulesh-comm.cpp
    ${lulesh_dir}/lulesh-init.cpp
    ${lulesh_dir}/lulesh-util.cpp
    ${lulesh_dir}/lulesh-viz.cpp
    )
set(lulesh_exe_source
    ${lulesh_dir}/lulesh.cpp
    )

## Ignore omp pragmas when openMP is not defined
set_source_files_properties(${lulesh_lib_sources} PROPERTIES COMPILE_FLAGS " -Wno-unknown-pragmas " )
set_source_files_properties(${lulesh_exe_source}  PROPERTIES COMPILE_FLAGS " -Wno-unknown-pragmas " )

if ( ENABLE_MPI )
    make_library( LIBRARY_NAME 
                    meshapiLuleshLib
                  LIBRARY_SOURCES 
                    "${lulesh_lib_sources}" 
                    "${lulesh_lib_headers}"
                  WITH_MPI
                  ) 
    
    add_dependencies(meshapiLuleshLib meshapi)
    target_link_libraries( meshapiLuleshLib  meshapi )
                      
    make_executable( EXECUTABLE_NAME 
                        meshapiLulesh 
                     EXECUTABLE_SOURCE 
                        ${lulesh_exe_source}
                     DEPENDS_ON 
                        "meshapi" 
                        "meshapiLuleshLib"
                     WITH_MPI
                   )
                       
else()
    make_library( LIBRARY_NAME 
                    meshapiLuleshLib
                  LIBRARY_SOURCES 
                    "${lulesh_lib_sources}" 
                    "${lulesh_lib_headers}"
                  ) 
    
    add_dependencies(meshapiLuleshLib meshapi)
    target_link_libraries( meshapiLuleshLib  meshapi )
                      
    make_executable( EXECUTABLE_NAME 
                        meshapiLulesh 
                     EXECUTABLE_SOURCE 
                        ${lulesh_exe_source}
                     DEPENDS_ON 
                        "meshapi" 
                        "meshapiLuleshLib"
                   )
endif()



# ---- Placed here temporarily --  determines whether to build the original lulesh 2.0 code
set(buildOrigLuleshMPI FALSE)
if ( ${buildOrigLuleshMPI} AND ENABLE_MPI)
    set(lulesh_dir lulesh2.0.3_orig  )

    set(lulesh_lib_headers
        ${lulesh_dir}/lulesh.h
        # ${lulesh_dir}/lulesh_tuple.hpp
        )

    set(lulesh_lib_sources
        ${lulesh_dir}/lulesh-comm.cc
        ${lulesh_dir}/lulesh-init.cc
        ${lulesh_dir}/lulesh-util.cc
        ${lulesh_dir}/lulesh-viz.cc
        )
    set(lulesh_exe_source
        ${lulesh_dir}/lulesh.cc
        )
    
    set_source_files_properties(${lulesh_lib_sources} PROPERTIES COMPILE_FLAGS " -Wno-unknown-pragmas -Wno-unused-parameter" )
    set_source_files_properties(${lulesh_exe_source}  PROPERTIES COMPILE_FLAGS " -Wno-unknown-pragmas -Wno-unused-parameter" )

    make_library( LIBRARY_NAME meshapiLuleshOrigLib LIBRARY_SOURCES "${lulesh_lib_sources}"  "${lulesh_lib_headers}" WITH_MPI )
    make_executable( EXECUTABLE_NAME meshapiOrigLulesh EXECUTABLE_SOURCE ${lulesh_exe_source} DEPENDS_ON "meshapiLuleshOrigLib" WITH_MPI )
endif()
                   