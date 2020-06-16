# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# 3rd Party Dependencies
#------------------------------------------------------------------------------

# Policy to use <PackageName>_ROOT variable in find_<Package> commands
# Policy added in 3.12+
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

#------------------------------------------------------------------------------
# UMPIRE
#------------------------------------------------------------------------------
if (UMPIRE_DIR)
    if (NOT EXISTS "${UMPIRE_DIR}")
        message(FATAL_ERROR "Given UMPIRE_DIR does not exist: ${UMPIRE_DIR}")
    endif()

    if (NOT IS_DIRECTORY "${UMPIRE_DIR}")
        message(FATAL_ERROR "Given UMPIRE_DIR is not a directory: ${UMPIRE_DIR}")
    endif()

    find_package(umpire REQUIRED PATHS ${UMPIRE_DIR} )

    message(STATUS "Checking for expected Umpire target 'umpire'")
    if (NOT TARGET umpire)
        message(FATAL_ERROR "Umpire failed to load: ${UMPIRE_DIR}")
    else()
        message(STATUS "Umpire loaded: ${UMPIRE_DIR}")
        set_property(TARGET umpire 
                     APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                    ${UMPIRE_INCLUDE_DIRS})
        set(UMPIRE_FOUND TRUE CACHE BOOL "")
    endif()
else()
    message(STATUS "Umpire support is OFF")
    set(UMPIRE_FOUND FALSE CACHE BOOL "")
endif()


#------------------------------------------------------------------------------
# RAJA
#------------------------------------------------------------------------------
if (RAJA_DIR)
    if (NOT EXISTS "${RAJA_DIR}")
        message(FATAL_ERROR "Given RAJA_DIR does not exist: ${RAJA_DIR}")
    endif()

    if (NOT IS_DIRECTORY "${RAJA_DIR}")
        message(FATAL_ERROR "Given RAJA_DIR is not a directory: ${RAJA_DIR}")
    endif()

    find_package(RAJA REQUIRED PATHS ${RAJA_DIR} )

    message(STATUS "Checking for expected RAJA target 'RAJA'")
    if (NOT TARGET RAJA)
        message(FATAL_ERROR "RAJA failed to load: ${RAJA_DIR}")
    else()
        message(STATUS "RAJA loaded: ${RAJA_DIR}")
        set(RAJA_FOUND TRUE CACHE BOOL "")
    endif()

    # Note: camp sets a compile feature that is not available on XL
    set_target_properties(camp PROPERTIES INTERFACE_COMPILE_FEATURES "")
else()
    message(STATUS "RAJA support is OFF" )
    set(RAJA_FOUND FALSE CACHE BOOL "")
endif()


#------------------------------------------------------------------------------
# HDF5
#------------------------------------------------------------------------------
if (HDF5_DIR)
    include(cmake/thirdparty/SetupHDF5.cmake)
else()
    message(STATUS "HDF5 support is OFF")
endif()


#------------------------------------------------------------------------------
# Conduit
#------------------------------------------------------------------------------
if (CONDUIT_DIR)
    include(cmake/thirdparty/FindConduit.cmake)

    # Manually set includes as system includes
    set_property(TARGET conduit::conduit 
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 "${CONDUIT_INSTALL_PREFIX}/include/")

    set_property(TARGET conduit::conduit 
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 "${CONDUIT_INSTALL_PREFIX}/include/conduit/")
else()
    message(STATUS "Conduit support is OFF")
endif()


#------------------------------------------------------------------------------
# MFEM
#------------------------------------------------------------------------------
if (MFEM_DIR)
    include(cmake/thirdparty/FindMFEM.cmake)
    blt_register_library( NAME      mfem
                          INCLUDES  ${MFEM_INCLUDE_DIRS}
                          LIBRARIES ${MFEM_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "MFEM support is OFF")
endif()


#------------------------------------------------------------------------------
# Shroud - Generates C/Fortran/Python bindings
#------------------------------------------------------------------------------
if(EXISTS ${SHROUD_EXECUTABLE})
    execute_process(COMMAND ${SHROUD_EXECUTABLE}
                    --cmake ${CMAKE_CURRENT_BINARY_DIR}/SetupShroud.cmake
                    ERROR_VARIABLE SHROUD_cmake_error
                    RESULT_VARIABLE SHROUD_cmake_result
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(NOT "${SHROUD_cmake_result}" STREQUAL "0")
        message(FATAL_ERROR "Error code from Shroud: ${SHROUD_cmake_result}\n${SHROUD_cmake_error}")
    endif()

    include(${CMAKE_CURRENT_BINARY_DIR}/SetupShroud.cmake)
else()
    message(STATUS "Shroud support is OFF")
endif()


#------------------------------------------------------------------------------
# SCR
#------------------------------------------------------------------------------
if (SCR_DIR)
    include(cmake/thirdparty/FindSCR.cmake)
    blt_register_library( NAME      scr
                          INCLUDES  ${SCR_INCLUDE_DIRS}
                          LIBRARIES ${SCR_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "SCR support is OFF")
endif()


# Remove exported OpenMP flags because they are not language agnostic
set(_props)
if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13.0" )
    list(APPEND _props INTERFACE_LINK_OPTIONS)
endif()
list(APPEND _props INTERFACE_COMPILE_OPTIONS)

foreach(_target RAJA camp umpire umpire_alloc)
    if(TARGET ${_target})
        message(STATUS "Removing OpenMP Flags from target[${_target}]")

        foreach(_prop ${_props})
            get_target_property(_flags ${_target} ${_prop})
            if ( _flags )
                string( REPLACE "${OpenMP_CXX_FLAGS}" ""
                        correct_flags "${_flags}" )
                string( REPLACE "${OpenMP_Fortran_FLAGS}" ""
                        correct_flags "${correct_flags}" )

                set_target_properties( ${_target} PROPERTIES ${_prop} "${correct_flags}" )
            endif()
        endforeach()
    endif()
endforeach()


#------------------------------------------------------------------------------
# LUA
#------------------------------------------------------------------------------
if (LUA_DIR)
    include(cmake/thirdparty/FindLUA.cmake)
    blt_register_library(
        NAME          lua
        INCLUDES      ${LUA_INCLUDE_DIR}
        LIBRARIES     ${LUA_LIBRARY}
        TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "LUA support is OFF")
    set(LUA_FOUND OFF CACHE BOOL "")
endif()

