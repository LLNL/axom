# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       Inlet
                        COMPONENTS SLIC Sidre )

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(inlet_headers
    Field.hpp
    Inlet.hpp
    Reader.hpp
    SchemaCreator.hpp
    Table.hpp 
    DocWriter.hpp
    SphinxDocWriter.hpp 
    inlet_utils.hpp )

set(inlet_sources
    Field.cpp
    Inlet.cpp
    Table.cpp
    SphinxDocWriter.cpp 
    inlet_utils.cpp )

blt_list_append(TO inlet_headers ELEMENTS LuaReader.hpp IF LUA_FOUND)
blt_list_append(TO inlet_sources ELEMENTS LuaReader.cpp IF LUA_FOUND)

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
set(inlet_depends core slic sidre fmt)

blt_list_append(TO inlet_depends ELEMENTS lua IF LUA_FOUND)

blt_add_library( NAME        inlet
                 SOURCES     ${inlet_sources}
                 HEADERS     ${inlet_headers}
                 DEPENDS_ON  ${inlet_depends}
                 FOLDER      axom/inlet
                 OBJECT      TRUE )

axom_write_unified_header(NAME    inlet
                          HEADERS ${inlet_headers} )

axom_install_component(NAME    inlet
                       HEADERS ${inlet_headers} )


#------------------------------------------------------------------------------
# Add examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(PREFIX inlet)
