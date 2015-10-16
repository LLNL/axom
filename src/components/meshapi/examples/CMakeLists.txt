
###########################
# Build the example sources
###########################
set(example_sources
	ShockTube.cpp
    UnstructMeshField.cpp
    TemplateEx.cpp
    )


#set_source_files_properties(TemplateEx.cpp PROPERTIES COMPILE_FLAGS "-S -fverbose-asm " )
#get_source_file_property(_origflags UnstructMeshField.cpp COMPILE_FLAGS)
#set_source_files_properties(UnstructMeshField.cpp "${_origflags} -fwhole-program -funroll-loops ")

foreach(example_source ${example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    make_executable(EXECUTABLE_NAME "meshapi${example_name}" EXECUTABLE_SOURCE ${example_source} DEPENDS_ON meshapi slic)
endforeach()


####################################
# Build meshapi version of lulesh 2.0
####################################

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
if ( NOT ENABLE_OMP )
    set_source_files_properties(${lulesh_lib_sources} ${lulesh_exe_source} PROPERTIES COMPILE_FLAGS "${DISABLE_OMP_PRAGMA_WARNINGS}" )
endif()


make_library( LIBRARY_NAME      meshapiLuleshLib
              LIBRARY_SOURCES   "${lulesh_lib_sources}" "${lulesh_lib_headers}"
              USE_OPENMP ${ENABLE_OMP} 
              USE_MPI ${ENABLE_MPI}              
              ) 

add_dependencies(meshapiLuleshLib meshapi)
target_link_libraries( meshapiLuleshLib  meshapi )
                  
make_executable( EXECUTABLE_NAME    meshapiLulesh 
                 EXECUTABLE_SOURCE  ${lulesh_exe_source}
                 DEPENDS_ON 
                    slic
                    meshapi 
                    meshapiLuleshLib
                  USE_OPENMP ${ENABLE_OMP} 
                  USE_MPI ${ENABLE_MPI}              
                  ) 


####################################
# Build original lulesh example
####################################

# ---- Option to control if original Lulesh builds
set(buildOrigLulesh TRUE)

if ( ${buildOrigLulesh} )
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
    
    # Add some extra configuration flags for MPI, OpenMP and allowing unused parameters (in 3rd party source code)
    set(luleshExtraCompileFlags ${DISABLE_UNUSED_PARAMETER_WARNINGS})
    if ( NOT ENABLE_OMP )
        set(luleshExtraCompileFlags "${luleshExtraCompileFlags} ${DISABLE_OMP_PRAGMA_WARNINGS}")
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
        LIBRARY_NAME meshapiLuleshOrigLib 
        LIBRARY_SOURCES "${lulesh_lib_sources}"  "${lulesh_lib_headers}"  
        USE_OPENMP ${ENABLE_OMP} 
        USE_MPI ${ENABLE_MPI}              
       ) 
    make_executable( 
        EXECUTABLE_NAME meshapiOrigLulesh 
        EXECUTABLE_SOURCE ${lulesh_exe_source} 
        DEPENDS_ON "meshapiLuleshOrigLib" 
        USE_OPENMP ${ENABLE_OMP} 
        USE_MPI ${ENABLE_MPI}              
       ) 
          
endif()



####################################
# Build meshapi version of tinyHydro
####################################
      
set(tinyHydro_dir tinyHydro  )

set(tinyHydro_lib_headers
    ${tinyHydro_dir}/TinyHydroTypes.hpp
    ${tinyHydro_dir}/HydroC.hpp
    ${tinyHydro_dir}/Part.hpp
    ${tinyHydro_dir}/State.hpp
    ${tinyHydro_dir}/VectorXY.hpp
    ${tinyHydro_dir}/PolygonMeshXY.hpp
    ${tinyHydro_dir}/myTimer.hpp        
    )

set(tinyHydro_lib_sources
    ${tinyHydro_dir}/TinyHydroTypes.cpp
    ${tinyHydro_dir}/HydroC.cpp
    ${tinyHydro_dir}/Part.cpp
    ${tinyHydro_dir}/State.cpp
    ${tinyHydro_dir}/PolygonMeshXY.cpp
    ${tinyHydro_dir}/myTimer.cpp        
    )


make_library( 
    LIBRARY_NAME meshapiTinyHydroLib 
    LIBRARY_SOURCES "${tinyHydro_lib_sources}"  "${tinyHydro_lib_headers}"
    DEPENDS_ON 
        slic
        meshapi 
   ) 
   
   
#
# Add gtest based tests
#
set(gtest_meshapi_tinyHydro_tests
    ${tinyHydro_dir}/tests/meshapi_tinyHydro_unitTests.cpp
    ${tinyHydro_dir}/tests/meshapi_tinyHydro_sedovTwoPart.cpp
    ${tinyHydro_dir}/tests/meshapi_tinyHydro_sod1DTwoPart.cpp
)


if(UNIX AND NOT APPLE)
    set(rtDep "-lrt")
else()
    set(rtDep "")
endif()
    
foreach(test ${gtest_meshapi_tinyHydro_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON meshapiTinyHydroLib slic ${rtDep})
endforeach()


# Add some extra flags to disable warnings about variables used only for SLIC_ASSERTS
# Temporarily used in inaccessible code block to get around warning.
# TODO: Replace with new macros from George..
#SET(EXTRA_WARNING_DISABLE_FLAGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable ")
#get_source_file_property(_origflags ${tinyHydro_dir}/tests/meshapi_tinyHydro_unitTests.cpp COMPILE_FLAGS)
#set_source_files_properties(
#         ${gtest_meshapi_tinyHydro_tests}
#         PROPERTIES COMPILE_FLAGS "${_origflags} ${EXTRA_WARNING_DISABLE_FLAGS}"     
#         )


