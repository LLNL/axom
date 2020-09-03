# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
axom_component_requires(NAME       Klee
                        TPLS       Conduit )

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(klee_headers
        Geometry.hpp
        IO.hpp
        Shape.hpp
        ShapeSet.hpp
    )

set(klee_sources
        Geometry.cpp
        IO.cpp
        Shape.cpp
        ShapeSet.cpp
    )

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
set(klee_depends core conduit::conduit)

blt_add_library(NAME        klee
                SOURCES     ${klee_sources}
                HEADERS     ${klee_headers}
                DEPENDS_ON  ${klee_depends}
                FOLDER      axom/klee
                OBJECT      TRUE
                )

axom_write_unified_header( NAME klee
                           HEADERS ${klee_headers}
                           )

axom_install_component( NAME    klee
                        HEADERS ${klee_headers}
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
# axom_add_code_checks(PREFIX klee)
