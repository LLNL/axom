if(ENABLE_PYTHON_API)
    set(ENABLE_SHARED_LIBS TRUE)
endif(ENABLE_PYTHON_API)

## add components

# Note:
# Shroud has to be enabled if we are building with fortran.  Some components need the shroud runtime library.
blt_add_component(COMPONENT_NAME shroud DEFAULT_STATE ON)

add_component(COMPONENT_NAME common DEFAULT_STATE ON)

# Lumberjack is an extension to SLIC that supports message filtering and reduction capabilities for parallel (MPI) executables.
# It is not meant for serial (non-mpi) executables.
if(ENABLE_MPI)
	blt_add_component(COMPONENT_NAME lumberjack DEFAULT_STATE ON)
endif(ENABLE_MPI)

blt_add_component(COMPONENT_NAME slic DEFAULT_STATE ON)
blt_add_component(COMPONENT_NAME slam DEFAULT_STATE ON)
blt_add_component(COMPONENT_NAME sidre DEFAULT_STATE ON)
blt_add_component(COMPONENT_NAME quest DEFAULT_STATE ON)

if(ENABLE_MPI)
    blt_add_component(COMPONENT_NAME spio DEFAULT_STATE ON)
endif(ENABLE_MPI)
