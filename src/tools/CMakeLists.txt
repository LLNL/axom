# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# This directory contains source code for utility executables that utilize
# Axom or that support Axom
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# convert_protocol is a utility that converts
# a datastore from the sidre_hdf5 protocol
# to any of the other supported protocols.
# It requires slam, sidre, mpi and hdf5.
#------------------------------------------------------------------------------
if( AXOM_ENABLE_SLAM AND AXOM_ENABLE_SIDRE AND ENABLE_MPI AND HDF5_FOUND )
    set(datastore_converter_sources
        convert_sidre_protocol.cpp
    )
    set(datastore_converter_depends 
        axom
        mpi
        cli11
        fmt)

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
# mesh_tester is a utility that tests a mesh for a variety of problems.
#------------------------------------------------------------------------------
set(mesh_tester_sources
    mesh_tester.cpp
)
set(mesh_tester_depends 
    axom
    )

blt_list_append(TO mesh_tester_depends ELEMENTS mpi IF ENABLE_MPI)
blt_list_append(TO mesh_tester_depends ELEMENTS openmp IF ENABLE_OPENMP)
blt_list_append(TO mesh_tester_depends ELEMENTS umpire IF UMPIRE_FOUND)
blt_list_append(TO mesh_tester_depends ELEMENTS RAJA IF RAJA_FOUND)
blt_list_append(TO mesh_tester_depends ELEMENTS cuda IF ENABLE_CUDA)

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
