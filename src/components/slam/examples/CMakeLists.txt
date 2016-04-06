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
    make_executable(NAME "slam_${example_name}" 
                    SOURCES ${example_source} 
                    DEPENDS_ON slam slic
                    ADD_CTEST IS_EXAMPLE
                    )
endforeach()


####################################
# Build slam (MeshAPI) version of lulesh 2.0
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
if ( NOT ENABLE_OPENMP )
    set_source_files_properties(${lulesh_lib_sources} ${lulesh_exe_source} 
        PROPERTIES COMPILE_FLAGS "${ATK_DISABLE_OMP_PRAGMA_WARNINGS}" )
endif()

make_library( NAME      slamLuleshLib
              SOURCES   "${lulesh_lib_sources}" "${lulesh_lib_headers}"
              USE_OPENMP ${ENABLE_OPENMP} 
              ) 

add_dependencies(slamLuleshLib slam)
target_link_libraries( slamLuleshLib  slam )
                  
make_executable( NAME    slam_Lulesh 
                 SOURCES  ${lulesh_exe_source}
                 DEPENDS_ON 
                    slic
                    slam 
                    slamLuleshLib
                  USE_OPENMP ${ENABLE_OPENMP} 
                  ADD_CTEST IS_EXAMPLE
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
        NAME slamLuleshOrigLib 
        SOURCES "${lulesh_lib_sources}"  "${lulesh_lib_headers}"  
        USE_OPENMP ${ENABLE_OPENMP} 
       ) 
    make_executable( 
        NAME slam_OrigLulesh 
        SOURCES ${lulesh_exe_source} 
        DEPENDS_ON "slamLuleshOrigLib" 
        USE_OPENMP ${ENABLE_OPENMP}
        IS_EXAMPLE
       ) 
          
endif()



####################################
# Build slam (MeshAPI) version of tinyHydro
####################################
      
set(tinyHydro_dir tinyHydro  )

set(tinyHydro_lib_headers
    ${tinyHydro_dir}/TinyHydroTypes.hpp
    ${tinyHydro_dir}/HydroC.hpp
    ${tinyHydro_dir}/Part.hpp
    ${tinyHydro_dir}/State.hpp
    ${tinyHydro_dir}/VectorXY.hpp
    ${tinyHydro_dir}/PolygonMeshXY.hpp
    )

set(tinyHydro_lib_sources
    ${tinyHydro_dir}/TinyHydroTypes.cpp
    ${tinyHydro_dir}/HydroC.cpp
    ${tinyHydro_dir}/Part.cpp
    ${tinyHydro_dir}/State.cpp
    ${tinyHydro_dir}/PolygonMeshXY.cpp
    )


make_library( 
    NAME slamTinyHydroLib 
    SOURCES "${tinyHydro_lib_sources}"  "${tinyHydro_lib_headers}"
    DEPENDS_ON 
        slic
        slam
   ) 
   
   
#
# Add gtest based tests
#
set(gtest_slam_tinyHydro_tests
    ${tinyHydro_dir}/tests/slam_tinyHydro_unitTests.cpp
)
foreach(test ${gtest_slam_tinyHydro_tests})
    add_gtest( TEST_SOURCE ${test} DEPENDS_ON slamTinyHydroLib slic)
endforeach()

make_executable(
    NAME "slam_tinyHydro_sedov" 
    SOURCES ${tinyHydro_dir}/tests/slam_tinyHydro_sedovTwoPart.cpp
    DEPENDS_ON slamTinyHydroLib
    IS_EXAMPLE
)
make_executable( 
    NAME "slam_tinyHydro_sod" 
    SOURCES ${tinyHydro_dir}/tests/slam_tinyHydro_sod1DTwoPart.cpp
    DEPENDS_ON slamTinyHydroLib
    IS_EXAMPLE
    ADD_CTEST
)
    




# Add some extra flags to disable warnings about variables used only for SLIC_ASSERTS
# Temporarily used in inaccessible code block to get around warning.
# TODO: Replace with new macros from George..
#SET(EXTRA_WARNING_DISABLE_FLAGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable ")
#get_source_file_property(_origflags ${tinyHydro_dir}/tests/slam_tinyHydro_unitTests.cpp COMPILE_FLAGS)
#set_source_files_properties(
#         ${gtest_slam_tinyHydro_tests}
#         PROPERTIES COMPILE_FLAGS "${_origflags} ${EXTRA_WARNING_DISABLE_FLAGS}"     
#         )


