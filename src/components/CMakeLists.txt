if(ENABLE_PYTHON_API)
    set(ENABLE_SHARED_LIBS TRUE)
endif(ENABLE_PYTHON_API)

## add components
## Note: Default state of components is controlled by the ENABLE_ALL_COMPONENTS
## option. However, the 'common' component is always enabled.

if ( ENABLE_SHROUD )
  blt_add_component(COMPONENT_NAME shroud DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
endif()

blt_add_component(COMPONENT_NAME axom_utils DEFAULT_STATE ON)

# Lumberjack is an extension to SLIC that supports message filtering and
# reduction capabilities for parallel (MPI) executables.It is not meant for
# serial (non-mpi) executables.
if(ENABLE_MPI)
  blt_add_component( COMPONENT_NAME lumberjack
                     DEFAULT_STATE  ${ENABLE_ALL_COMPONENTS})
endif()

blt_add_component(COMPONENT_NAME slic DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
blt_add_component(COMPONENT_NAME primal DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
blt_add_component(COMPONENT_NAME slam DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
blt_add_component(COMPONENT_NAME sidre DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
blt_add_component(COMPONENT_NAME mint DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})

if ( BOOST_FOUND )
  blt_add_component(COMPONENT_NAME quest DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
endif()

if(ENABLE_MPI)
  blt_add_component(COMPONENT_NAME spio DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
endif()


## Add dependency for each enabled component on the copy_axom_includes target
## to ensure that the include files are copied over
## Note: This assumes components have a cmake target with the same name 
## (e.g. that slic has a target named 'slic')
set(axom_components shroud axom_utils lumberjack slic primal slam sidre mint quest spio)
foreach(comp ${axom_components})
   if(NOT ${comp} STREQUAL shroud)  # temporary so don't indent
   string(TOUPPER ${comp} COMP)
   if(ENABLE_${COMP})
       add_dependencies(${comp} copy_axom_includes)
   endif()
   endif()
endforeach()
