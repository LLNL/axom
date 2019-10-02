# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

################################
# Check necessary dependencies
################################
axom_component_requires(NAME       Curves
                        COMPONENTS SLIC)

#####################################
# Specify component headers/sources
#####################################
set( curves_headers
     curves.hpp
     )

set( curves_sources
     curves.cpp
     )


#####################################
# Specify dependencies
#####################################
set( curves_dependencies
     core
     slic
   )

## add any other dependencies here
## blt_list_append( TO curves_dependencies ELEMENTS raja IF ${RAJA_FOUND})
## etc.


################################
# Make/Install the library
################################
blt_add_library(
    NAME       curves
    SOURCES    ${curves_sources}
    HEADERS    ${curves_headers}
    DEPENDS_ON ${curves_dependencies}
    FOLDER     axom/curves
    OBJECT     TRUE
    )


axom_write_unified_header( NAME curves HEADERS ${curves_headers} )

axom_install_component( NAME curves HEADERS ${curves_headers} )

################################
# Add examples
################################
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX slic)