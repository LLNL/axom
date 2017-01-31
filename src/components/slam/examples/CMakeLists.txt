###############
# Slam examples
###############

#set_source_files_properties(TemplateEx.cpp PROPERTIES COMPILE_FLAGS "-S -fverbose-asm " )
#get_source_file_property(_origflags UnstructMeshField.cpp COMPILE_FLAGS)
#set_source_files_properties(UnstructMeshField.cpp "${_origflags} -fwhole-program -funroll-loops ")



#
# Single source standalone examples
#
set(example_sources
    ShockTube.cpp
    TemplateEx.cpp
    )

set(slam_tests_depends_on common slic slam)

# Slam works without boost, but we need to add the header dependency when boost is available
set(slam_tests_extra_depends)
if(ENABLE_BOOST)
  list(APPEND slam_tests_extra_depends boost)
endif()

foreach(example_source ${example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME        "slam_${example_name}"
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${slam_tests_depends_on} ${slam_tests_extra_depends}
        )
        
    blt_add_test(
        NAME        slam_${example_name} 
        COMMAND     slam_${example_name}
        )
endforeach()

#
# Slam (MeshAPI) version of lulesh 2.0
#
add_subdirectory(lulesh2.0.3)

# Original lulesh example

set(buildOrigLulesh TRUE)        # ---- Controls whether original Lulesh should build
if ( ${buildOrigLulesh} )
   add_subdirectory(lulesh2.0.3_orig)
endif()


#
# Slam (MeshAPI) version of tinyHydro
#
add_subdirectory(tinyHydro)


#
# An example that use the Slam iterators and therefore depends on boost
#
if(ENABLE_BOOST)
    list(APPEND slam_tests_extra_depends fmt)

    # Unstructured hexahedral mesh
    blt_add_executable(
        NAME        slam_unstructMesh
        SOURCES     UnstructMeshField.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${slam_tests_depends_on} ${slam_tests_extra_depends}
        )
        
    blt_add_test(
        NAME        slam_unstructMesh
        COMMAND     slam_unstructMesh ${CMAKE_CURRENT_SOURCE_DIR}/../data
        )

endif()
