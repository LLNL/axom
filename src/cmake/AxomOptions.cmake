# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Defines CMake options for Axom's build system
#------------------------------------------------------------------------------

option(AXOM_ENABLE_ANNOTATIONS "Enables code annotations to facilitate performance evaluation." OFF)
option(AXOM_ENABLE_SPARSEHASH "Enables Sparsehash." ON)
option(AXOM_ENABLE_ALL_COMPONENTS "Enables all components by default" ON)
option(AXOM_USE_64BIT_INDEXTYPE "Use 64-bit integers for axom::IndexType" OFF)
option(AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION "Enable Axom's version of the MFEM SidreDataCollection" OFF)

if(NOT CMAKE_CONFIGURATION_TYPES)
    if(CMAKE_BUILD_TYPE MATCHES "(Debug|RelWithDebInfo)")
        option(AXOM_ENABLE_EXPORTS "Add in symbols to demangle axom function names in stacktraces" ON)
    else()
    	option(AXOM_ENABLE_EXPORTS "Add in symbols to demangle axom function names in stacktraces" OFF)
    endif()
endif()

cmake_dependent_option(AXOM_ENABLE_TESTS "Enables Axom Tests" ON "ENABLE_TESTS" OFF)
cmake_dependent_option(AXOM_ENABLE_DOCS "Enables Axom Docs" ON "ENABLE_DOCS" OFF)
cmake_dependent_option(AXOM_ENABLE_EXAMPLES "Enables Axom Examples" ON "ENABLE_EXAMPLES" OFF)
option(AXOM_ENABLE_TOOLS "Enables Axom Tools" ON)

cmake_dependent_option(AXOM_ENABLE_MPI3 "Enables use of MPI-3 features" OFF "ENABLE_MPI" OFF)
mark_as_advanced(AXOM_ENABLE_MPI3)
