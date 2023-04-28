# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Add Axom components
#------------------------------------------------------------------------------
# Note: Default state of components is controlled by the
# AXOM_ENABLE_ALL_COMPONENTS option.
#------------------------------------------------------------------------------

# Core is an essential part of Axom and cannot be turned off
axom_add_component(COMPONENT_NAME core DEFAULT_STATE ON)
if(NOT AXOM_ENABLE_CORE)
    message(FATAL_ERROR "Core is an essential part of Axom and cannot be turned off")
endif()

# Lumberjack is a parallel message filtering and reduction library. It can be used
# by itself or as a SLIC Log Stream. It is not meant for serial (non-mpi) executables.
if(AXOM_ENABLE_MPI)
    axom_add_component(COMPONENT_NAME lumberjack
                       DEFAULT_STATE  ${AXOM_ENABLE_ALL_COMPONENTS})
else()
    message(STATUS "Axom Component Lumberjack turned off due to AXOM_ENABLE_MPI=OFF")
    set(AXOM_ENABLE_LUMBERJACK OFF CACHE BOOL "")
endif()

# Add components so that later ones depend on earlier ones
# (i.e. quest depends on mint, so quest follows mint)
axom_add_component(COMPONENT_NAME slic   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME slam   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME primal DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME sidre  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME mint   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME spin   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME inlet  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME klee   DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME quest  DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})
axom_add_component(COMPONENT_NAME multimat DEFAULT_STATE ${AXOM_ENABLE_ALL_COMPONENTS})

install(TARGETS              ${AXOM_COMPONENTS_ENABLED}
        EXPORT               axom-targets
        DESTINATION          lib)
install(EXPORT               axom-targets 
        NAMESPACE            axom::
        DESTINATION          lib/cmake)

#------------------------------------------------------------------------------
# Generate export symbols
#------------------------------------------------------------------------------
include(GenerateExportHeader)
set(_export_header ${PROJECT_BINARY_DIR}/axom_export_symbols)
generate_export_header(core
                       BASE_NAME axom
                       EXPORT_FILE_NAME ${_export_header})

#------------------------------------------------------------------------------
# Output some information about the configuration
#------------------------------------------------------------------------------
foreach(_component ${AXOM_COMPONENTS_FULL})
    string(TOUPPER ${_component} COMPONENT_NAME_UPPERCASE)
    if(AXOM_ENABLE_${COMPONENT_NAME_UPPERCASE})
        message(STATUS "Axom component ${_component} is ON")
    else()
        message(STATUS "Axom component ${_component} is OFF")
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
    set(_code_check_targets "docs" "doxygen_docs"
                            "check" "style" 
                            "clangformat_check" "clangformat_style" )
    foreach(_tgt ${_code_check_targets})
        if(TARGET ${_tgt})
            set_target_properties(${_tgt} PROPERTIES FOLDER "axom/code_checks")
        endif()
    endforeach()
endif()


#------------------------------------------------------------------------------
# Bugfix for dllimport/dllexport functionality in msvc with object libraries
# Each axom component needs to have 'axom_EXPORTS' defined when building
#------------------------------------------------------------------------------
if(WIN32)
    foreach(_component ${AXOM_COMPONENTS_ENABLED})
        blt_add_target_definitions(
            TO                 ${_component}
            SCOPE              PRIVATE
            TARGET_DEFINITIONS "axom_EXPORTS")
    endforeach()
endif()
