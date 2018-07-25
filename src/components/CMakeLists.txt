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

## add components
## Note: Default state of components is controlled by the
## AXOM_ENABLE_ALL_COMPONENTS option. However, the 'axom_utils' component is always
## enabled.

axom_add_component(COMPONENT_NAME axom_utils DEFAULT_STATE ON)

# Lumberjack is an extension to SLIC that supports message filtering and
# reduction capabilities for parallel (MPI) executables. It is not meant for
# serial (non-mpi) executables.
if(ENABLE_MPI)
    axom_add_component( COMPONENT_NAME lumberjack
                        DEFAULT_STATE  ${AXOM_ENABLE_ALL_COMPONENTS})
endif()

axom_add_component(COMPONENT_NAME slic DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME primal DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME slam DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME sidre DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME mint DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME quest DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})


## Add dependency for each enabled component on the copy_axom_includes target
## to ensure that the include files are copied over
## Note: This assumes components have a cmake target with the same name 
## (e.g. that slic has a target named 'slic')
set(axom_components axom_utils lumberjack slic primal slam sidre mint quest)
foreach(comp ${axom_components})
    string(TOUPPER ${comp} COMP)
    if(ENABLE_${COMP})
        add_dependencies(${comp} copy_axom_includes)
    endif()
endforeach()
