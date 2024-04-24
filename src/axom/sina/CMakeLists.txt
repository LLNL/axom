# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina -- API to write data to Sina's common file format
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       Sina
                        COMPONENTS SLIC
                        TPLS       Conduit )

#------------------------------------------------------------------------------
# Add the sina sources
#------------------------------------------------------------------------------
set(sina_headers
    include/sina/ConduitUtil.hpp
    include/sina/CppBridge.hpp
    include/sina/Curve.hpp
    include/sina/CurveSet.hpp
    include/sina/DataHolder.hpp
    include/sina/Datum.hpp
    include/sina/Document.hpp
    include/sina/File.hpp
    include/sina/ID.hpp
    include/sina/Record.hpp
    include/sina/Relationship.hpp
    include/sina/Run.hpp
    include/sina/sina.hpp )

set(sina_sources
    src/ConduitUtil.cpp
    src/Curve.cpp
    src/CurveSet.cpp
    src/DataHolder.cpp
    src/Datum.cpp
    src/Document.cpp
    src/File.cpp
    src/ID.cpp
    src/Record.cpp
    src/Relationship.cpp
    src/Run.cpp )

# # Add Adiak header and source
# if (AXOM_BUILD_SINA_ADIAK_BINDINGS)  # BRIAN LOOK AT ME RIGHT HERE -> you modified this variable at the end of the week. Is this correct? Should it be modified in the src code like you did?
#     list(APPEND sina_headers include/sina/AdiakWriter.hpp)
#     list(APPEND sina_sources src/AdiakWriter.cpp)
# endif()

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(sina_depends
    slic
    conduit::conduit
    )

# if(AXOM_BUILD_SINA_ADIAK_BINDINGS)
#     blt_list_append(TO sidre_depends ELEMENTS adiak)
# endif()

axom_add_library(NAME       sina
                 SOURCES    ${sina_sources}
                 HEADERS    ${sina_headers}
                 DEPENDS_ON ${sina_depends}
                 FOLDER     axom/sina)

axom_write_unified_header(NAME    sina
                          HEADERS ${sina_headers})

axom_install_component(NAME    sina
                       HEADERS ${sina_headers})

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if(AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
