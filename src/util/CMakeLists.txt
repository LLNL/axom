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

#------------------------------------------------------------------------------
# This directory contains source code for utility executables that utilize
# Axom or that support Axom
#------------------------------------------------------------------------------


##################################################
# convert_protocol is a utility that converts
# a datastore from the sidre_hdf5 protocol
# to any of the other supported protocols.
# It requires slam, sidre, mpi and hdf5.
##################################################
if( ENABLE_SLAM AND ENABLE_SIDRE AND ENABLE_MPI AND HDF5_FOUND )
    set(datastore_converter_sources
        convert_sidre_protocol.cpp
    )
    set(datastore_converter_depends 
        sidre conduit conduit_relay hdf5
        slic fmt 
        slam 
        mpi
        )

    if(ENABLE_LUMBERJACK)
      list(APPEND datastore_converter_depends lumberjack)
    endif()


    if(BOOST_FOUND)
      list(APPEND datastore_converter_depends boost)
    endif()

    if(SCR_FOUND)
      list(APPEND datastore_converter_depends scr)
    endif()

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

##################################################
# mesh_tester is a utility that tests a mesh for
# a variety of problems.
##################################################
set(mesh_tester_sources
    mesh_tester.cpp
)
set(mesh_tester_depends 
    axom_utils slic mint slam quest primal
    )


if(ENABLE_QUEST)
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

################################
# Add code checks
################################
axom_add_code_checks(PREFIX util)
