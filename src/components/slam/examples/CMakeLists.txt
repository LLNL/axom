#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

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
    PolicyPrototype.cpp
    )

set(slam_tests_depends_on slic slam)

set(slam_tests_extra_depends)

foreach(example_source ${example_sources})
    get_filename_component(example_name ${example_source} NAME_WE)
    blt_add_executable(
        NAME        "slam_${example_name}_ex"
        SOURCES     ${example_source}
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${slam_tests_depends_on} ${slam_tests_extra_depends}
        FOLDER      axom/slam/examples
        )

    blt_add_test(
        NAME        slam_${example_name}
        COMMAND     slam_${example_name}_ex
        )
endforeach()


#
# Slam (MeshAPI) version of tinyHydro
#
add_subdirectory(tinyHydro)


#
# An example that use the Slam iterators and therefore depends on C++11
#
if( NOT BLT_CXX_STD STREQUAL c++98 ) 
    list(APPEND slam_tests_extra_depends fmt)

    # Unstructured hexahedral mesh
    blt_add_executable(
        NAME        slam_unstructMesh_ex
        SOURCES     UnstructMeshField.cpp
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${slam_tests_depends_on} ${slam_tests_extra_depends}
        FOLDER      axom/slam/examples
        )

    blt_add_test(
        NAME        slam_unstructMesh
        COMMAND     slam_unstructMesh_ex ${CMAKE_CURRENT_SOURCE_DIR}/../data
        )

endif()


#
# Slam (MeshAPI) version of lulesh 2.0
# Note: Disabled on windows builds due to usage of unix APIs.
#
if(NOT WIN32)
	add_subdirectory(lulesh2.0.3)

	# Original lulesh example

	set(buildOrigLulesh TRUE)        # ---- Controls whether original Lulesh should build
	if ( ${buildOrigLulesh} )
	   add_subdirectory(lulesh2.0.3_orig)
	endif()
endif()
