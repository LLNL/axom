# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of examples
#------------------------------------------------------------------------------

if (LUA_FOUND)
    blt_add_executable( NAME       inlet_array_example 
                        SOURCES    arrays.cpp
                        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet
                        FOLDER     axom/inlet/examples )

    blt_add_executable( NAME       inlet_verification_example 
                        SOURCES    verification.cpp
                        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet
                        FOLDER     axom/inlet/examples )
    blt_add_test( NAME    inlet_verification 
                  COMMAND inlet_verification_example )

    blt_add_executable( NAME       inlet_documentation_generation_example 
                        SOURCES    documentation_generation.cpp
                        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                        DEPENDS_ON inlet
                        FOLDER     axom/inlet/examples )
    blt_add_test( NAME    inlet_documentation_generation 
                  COMMAND inlet_documentation_generation_example --enableDocs --file ${CMAKE_CURRENT_LIST_DIR}/example1.lua )

endif()
