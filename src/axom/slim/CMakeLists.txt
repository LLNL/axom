# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# SLIM component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(slim_headers
    Field.hpp
    GroupField.hpp
    IntField.hpp
    Map.hpp
    Structure.hpp
    )

set(slim_sources
    Structure.cpp
    )


blt_list_append(TO slim_headers ELEMENTS LuaMap.hpp IF LUA_FOUND)
blt_list_append(TO slim_sources ELEMENTS LuaMap.cpp IF LUA_FOUND)

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
set(slim_depends core slic fmt)

blt_list_append(TO slim_depends ELEMENTS lua IF LUA_FOUND)

blt_add_library( NAME        slim
                 SOURCES     ${slim_sources}
                 HEADERS     ${slim_headers}
                 DEPENDS_ON  ${slim_depends}
                 FOLDER      axom/slim
                 OBJECT      TRUE
                 )

axom_write_unified_header(NAME    slim
                          HEADERS ${slim_headers})

axom_install_component(NAME      slim
                       HEADERS   ${slim_headers}
                       )


#------------------------------------------------------------------------------
# Add examples
#------------------------------------------------------------------------------
#if (AXOM_ENABLE_EXAMPLES)
#  add_subdirectory(examples)
#endif()

#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(PREFIX slim)
