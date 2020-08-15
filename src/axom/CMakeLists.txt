# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Add Axom components
#------------------------------------------------------------------------------
# Note: Default state of components is controlled by the
# AXOM_ENABLE_ALL_COMPONENTS option.
#------------------------------------------------------------------------------

# Core is an essential part of Axom and cannot be turned off
add_subdirectory(core)

# Lumberjack is a parallel message filtering and reduction library. It can be used
# by itself or as a SLIC Log Stream. It is not meant for serial (non-mpi) executables.
if(ENABLE_MPI)
    axom_add_component( COMPONENT_NAME lumberjack
                        DEFAULT_STATE  ${AXOM_ENABLE_ALL_COMPONENTS})
else()
    message(STATUS "Axom Component Lumberjack turned off due to ENABLE_MPI=OFF")
    set(AXOM_ENABLE_LUMBERJACK OFF CACHE BOOL "")
endif()

# Add components so that later ones depend on earlier ones
# (i.e. quest depends on mint, so quest follows mint)
axom_add_component(COMPONENT_NAME slic   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME slam   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME mint   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME primal DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME spin   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME sidre  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
if (AXOM_ENABLE_SIDRE)
    axom_add_component(COMPONENT_NAME inlet  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
else()
    message(STATUS "Axom Component Inlet turned off due to AXOM_ENABLE_SIDRE=OFF")
    set(AXOM_ENABLE_INLET OFF CACHE BOOL "")
endif()
axom_add_component(COMPONENT_NAME quest  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})

# Combine all component object libraries into a unified library
blt_add_library(NAME       axom
                SOURCES    Axom.cpp
                DEPENDS_ON ${AXOM_COMPONENTS_ENABLED})

install(TARGETS              axom
        EXPORT               axom-targets
        DESTINATION          lib)
install(EXPORT axom-targets DESTINATION lib/cmake)


#------------------------------------------------------------------------------
# Output some information about the configuration
#------------------------------------------------------------------------------
foreach(comp ${AXOM_COMPONENTS_FULL})
    string(TOUPPER ${comp} COMPONENT_NAME_UPPERCASE)
    if(AXOM_ENABLE_${COMPONENT_NAME_UPPERCASE})
        message(STATUS "Axom component ${comp} is ON")
    else()
        message(STATUS "Axom component ${comp} is OFF")
    endif()
endforeach()

if(AXOM_DATA_DIR)
    message(STATUS "AXOM_DATA_DIR: ${AXOM_DATA_DIR}")
else()
    message(STATUS "AXOM_DATA_DIR: <undefined>")
endif()

#------------------------------------------------------------------------------
# Fix FOLDER property for some targets
#------------------------------------------------------------------------------
if(TARGET axom)
    set_target_properties(axom PROPERTIES FOLDER "axom")
endif()

# If Axom is the main project, set FOLDER property for top-level check targets
if ("${PROJECT_SOURCE_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
    #set(_code_check_targets "docs" "doxygen_docs"
    set(_code_check_targets "docs"
                            "check" "style" 
                            "uncrustify_check" "uncrustify_style" )
    foreach(_tgt ${_code_check_targets})
        if(TARGET ${_tgt})
            set_target_properties(${_tgt} PROPERTIES FOLDER "axom/code_checks")
        endif()
    endforeach()
endif()
