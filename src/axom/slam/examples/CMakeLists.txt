# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Slam examples
#------------------------------------------------------------------------------
# The following flags might be helpful for analyzing performance:
#   for TemplateEx:        -S -fverbose-asm
#   for UnstructMeshField: -fwhole-program -funroll-loops
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Single source standalone examples
#------------------------------------------------------------------------------
set(example_sources
    ShockTube.cpp
    PolicyPrototype.cpp
    UserDocs.cpp
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
        FOLDER      axom/slam/examples )

    if(AXOM_ENABLE_TESTS)
        blt_add_test(
            NAME        slam_${example_name}
            COMMAND     slam_${example_name}_ex )
    endif()
endforeach()


#------------------------------------------------------------------------------
# Slam version of tinyHydro
#------------------------------------------------------------------------------
add_subdirectory(tinyHydro)


#------------------------------------------------------------------------------
# Unstructured hexahedral mesh example
#------------------------------------------------------------------------------
list(APPEND slam_tests_extra_depends fmt)

blt_add_executable(
    NAME        slam_unstructMesh_ex
    SOURCES     UnstructMeshField.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${slam_tests_depends_on} ${slam_tests_extra_depends}
    FOLDER      axom/slam/examples )

if(AXOM_ENABLE_TESTS AND AXOM_DATA_DIR)
    blt_add_test(
        NAME        slam_unstructMesh
        COMMAND     slam_unstructMesh_ex ${AXOM_DATA_DIR}/slam
        )
endif()


#------------------------------------------------------------------------------
# Slam version of lulesh 2.0
# Note: Disabled on windows builds due to Lulesh's usage of unix APIs.
#------------------------------------------------------------------------------
if(NOT WIN32)
    add_subdirectory(lulesh2.0.3)

    # Original lulesh example
    set(_buildOrigLulesh TRUE)          # ---- Should we build original Lulesh?
    if ( ${_buildOrigLulesh} )
       add_subdirectory(lulesh2.0.3_orig)
    endif()
    unset(_buildOrigLulesh)
endif()

