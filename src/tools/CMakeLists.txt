# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# This directory contains source code for utility executables that utilize
# Axom or that support Axom
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# convert_protocol is a utility that converts a datastore from the sidre_hdf5 
# protocol to any of the other supported protocols.
# It requires slam, sidre, mpi and hdf5.
#------------------------------------------------------------------------------
if( AXOM_ENABLE_SLAM AND AXOM_ENABLE_SIDRE AND ENABLE_MPI AND HDF5_FOUND )
    set(datastore_converter_sources
        convert_sidre_protocol.cpp )

    set(datastore_converter_depends 
        cli11
        sidre conduit::conduit hdf5
        slic fmt 
        slam 
        mpi)

    blt_list_append(TO datastore_converter_depends ELEMENTS lumberjack mpi IF AXOM_ENABLE_LUMBERJACK)
    blt_list_append(TO datastore_converter_depends ELEMENTS scr IF SCR_FOUND)

    blt_add_executable(
        NAME       convert_sidre_protocol
        SOURCES    ${datastore_converter_sources}
        OUTPUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON ${datastore_converter_depends}
        FOLDER     axom/examples
        )
    install(TARGETS              convert_sidre_protocol
            DESTINATION          bin
            )
endif()

#------------------------------------------------------------------------------
# verify_mesh_blueprint is a utility that reads in a sidre datastore
# and verifies that the supplied path coresponds to a valid mesh blueprint.
# It requires slam, sidre, mpi and hdf5.
#------------------------------------------------------------------------------
if( AXOM_ENABLE_SLAM AND AXOM_ENABLE_SIDRE AND ENABLE_MPI AND HDF5_FOUND )
    set(verify_blueprint_sources
        verify_mesh_blueprint.cpp )

    set(verify_blueprint_depends 
        cli11
        sidre conduit::conduit hdf5
        slic fmt 
        slam 
        mpi)

    blt_list_append(TO verify_blueprint_depends ELEMENTS lumberjack mpi IF AXOM_ENABLE_LUMBERJACK)
    blt_list_append(TO verify_blueprint_depends ELEMENTS scr IF SCR_FOUND)

    blt_add_executable(
        NAME       verify_mesh_blueprint
        SOURCES    ${verify_blueprint_sources}
        OUTPUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON ${verify_blueprint_depends}
        FOLDER     axom/examples
        )
    install(TARGETS              verify_mesh_blueprint
            DESTINATION          bin
            )
endif()

#------------------------------------------------------------------------------
# mesh_tester is a utility that tests a mesh for a variety of problems.
#------------------------------------------------------------------------------
set(mesh_tester_sources
    mesh_tester.cpp
)
set(mesh_tester_depends 
    core slic mint slam quest primal cli11
    )

blt_list_append(TO mesh_tester_depends ELEMENTS lumberjack mpi IF AXOM_ENABLE_LUMBERJACK)
blt_list_append(TO mesh_tester_depends ELEMENTS openmp IF ENABLE_OPENMP)

if(AXOM_ENABLE_QUEST)
    blt_add_executable(
        NAME        mesh_tester
        SOURCES     ${mesh_tester_sources}
        OUTPUT_DIR  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${mesh_tester_depends}
        FOLDER      axom/examples
        )
    install(TARGETS              mesh_tester
            DESTINATION          bin
            )
endif()

#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(PREFIX tools)
