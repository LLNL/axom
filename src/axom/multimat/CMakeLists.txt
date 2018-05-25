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

##################################################################################
# MultiMat
##################################################################################

################################
# Check necessary dependencies
################################
axom_component_requires(NAME       MultiMat
                        COMPONENTS SLAM SLIC)

################################
# Specify all headers/sources
################################
set(multimat_headers
    multimat.hpp
    )

set(multimat_sources
    multimat.cpp
    )

################################
# Make/Install the library
################################
set(multimat_depends_on core slam)

blt_add_library(
    NAME        multimat
    SOURCES     ${multimat_sources}
    HEADERS     ${multimat_headers}
    DEPENDS_ON  ${multimat_depends_on}
    FOLDER      axom/multimat
    )

axom_install_component(NAME    multimat
                       HEADERS ${multimat_headers}
                       )

################################
# Add examples
################################
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests and benchmarks
################################
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add code checks
################################
axom_add_code_checks(
        PREFIX multimat
        )
