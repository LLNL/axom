# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# MultiMat component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       MultiMat
                        COMPONENTS SLAM SLIC)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(multimat_headers
    multimat.hpp
    mmfield.hpp
    mmsubfield.hpp
    )

set(multimat_sources
    multimat.cpp
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
axom_add_library(
    NAME        multimat
    SOURCES     ${multimat_sources}
    HEADERS     ${multimat_headers}
    DEPENDS_ON  slam
    FOLDER      axom/multimat
    )

# Set file back to C++ due to nvcc compiler error
set_source_files_properties(multimat.cpp PROPERTIES LANGUAGE CXX)

axom_write_unified_header(NAME    multimat
                          HEADERS ${multimat_headers}
                          )

axom_install_component(NAME    multimat
                       HEADERS ${multimat_headers}
                       )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()


#------------------------------------------------------------------------------
# Add docs and code checks
#------------------------------------------------------------------------------
axom_add_code_checks(
        PREFIX multimat
        )
