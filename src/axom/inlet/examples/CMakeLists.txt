# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Tests for Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

# Add Serial GTests based tests
if (LUA_FOUND)
    set(inlet_examples
        document_generation_example.cpp )

    foreach(ex ${inlet_examples})
        get_filename_component( ex_name ${ex} NAME_WE )
        blt_add_executable( NAME       ${ex_name}_example
                            SOURCES    ${ex}
                            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                            DEPENDS_ON inlet
                            FOLDER     axom/inlet/examples )
        blt_add_test( NAME    ${ex_name} 
                      COMMAND ${ex_name}_example )
    endforeach()
endif()
