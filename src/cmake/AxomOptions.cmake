#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

option(AXOM_ENABLE_SPARSEHASH "Enables Sparsehash." ON)
option(AXOM_ENABLE_ALL_COMPONENTS "Enables all components by default" ON)

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

cmake_dependent_option(AXOM_USE_MPI3 "Enables use of MPI-3 features" OFF "ENABLE_MPI" OFF)
mark_as_advanced(AXOM_USE_MPI3)
