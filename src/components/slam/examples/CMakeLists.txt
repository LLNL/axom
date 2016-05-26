###############
# Slam examples
###############

#
# Single source standalone examples
#
set(example_sources
        ShockTube.cpp
        UnstructMeshField.cpp
        TemplateEx.cpp
        )


#set_source_files_properties(TemplateEx.cpp PROPERTIES COMPILE_FLAGS "-S -fverbose-asm " )
#get_source_file_property(_origflags UnstructMeshField.cpp COMPILE_FLAGS)
#set_source_files_properties(UnstructMeshField.cpp "${_origflags} -fwhole-program -funroll-loops ")


# blt_add_executable depends_on handles our include dir settings
# these teste are using boost so boost needs to be a dep
set(slam_tests_depends_on common slic slam)
if (ENABLE_BOOST)
    list(APPEND slam_tests_depends_on boost)
endif()


foreach(example_source ${example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME        "slam_${example_name}"
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${slam_tests_depends_on}
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



#
# Original lulesh example
#

# ---- Option to control if original Lulesh builds
set(buildOrigLulesh TRUE)

if ( ${buildOrigLulesh} )
   add_subdirectory(lulesh2.0.3_orig)
endif()



#
# Slam (MeshAPI) version of tinyHydro
#
add_subdirectory(tinyHydro)
