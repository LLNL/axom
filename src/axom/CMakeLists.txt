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

## Add components
## Note: Default state of components is controlled by the
## AXOM_ENABLE_ALL_COMPONENTS option.

# Core is an essential part of Axom and cannot be turned off
add_subdirectory(core)

# Lumberjack is a parallel message filtering and reduction library. It can be used
# by itself or as a SLIC Log Stream. It is not meant for serial (non-mpi) executables.
if(ENABLE_MPI)
    axom_add_component( COMPONENT_NAME lumberjack
                        DEFAULT_STATE  ${AXOM_ENABLE_ALL_COMPONENTS})
else()
    set(AXOM_ENABLE_LUMBERJACK OFF CACHE BOOL "")
endif()

axom_add_component(COMPONENT_NAME slic   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME primal DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME slam   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME sidre  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME mint   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME quest  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})

foreach(comp ${AXOM_COMPONENTS_FULL})
    string(TOUPPER ${comp} COMPONENT_NAME_UPPERCASE)
    if(AXOM_ENABLE_${COMPONENT_NAME_UPPERCASE})
        message(STATUS "Axom component ${comp} is ON")
    else()
        message(STATUS "Axom component ${comp} is OFF")
    endif()
endforeach()
