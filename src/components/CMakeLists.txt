if(ENABLE_PYTHON_MODULE)
    set(BUILD_SHARED_LIBS TRUE)
endif(ENABLE_PYTHON_MODULE)

add_subdirectory(shroud)

## add components
add_component(COMPONENT_NAME common DEFAULT_STATE ON)
add_component(COMPONENT_NAME slic DEFAULT_STATE ON)
add_component(COMPONENT_NAME meshapi DEFAULT_STATE ON)
add_component(COMPONENT_NAME sidre DEFAULT_STATE ON)


if(ENABLE_MPI)
	add_component(COMPONENT_NAME lumberjack DEFAULT_STATE ON)
endif(ENABLE_MPI)
