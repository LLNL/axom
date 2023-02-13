# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
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

set(TPL_DEPS)

#------------------------------------------------------------------------------
# Create global variable to toggle between GPU targets
#------------------------------------------------------------------------------
if(AXOM_ENABLE_CUDA)
    set(axom_device_depends cuda CACHE STRING "" FORCE)
endif()
if(AXOM_ENABLE_HIP)
    set(axom_device_depends blt::hip CACHE STRING "" FORCE)
endif()

#------------------------------------------------------------------------------
# Camp (needed by RAJA and Umpire)
#------------------------------------------------------------------------------
if ((RAJA_DIR OR UMPIRE_DIR) AND NOT CAMP_DIR)
    message(FATAL_ERROR "CAMP_DIR is required if RAJA_DIR or UMPIRE_DIR is provided.")
endif()

if (CAMP_DIR)
    if (NOT EXISTS "${CAMP_DIR}")
        message(FATAL_ERROR "Given CAMP_DIR does not exist: ${CAMP_DIR}")
    endif()

    if (NOT IS_DIRECTORY "${CAMP_DIR}")
        message(FATAL_ERROR "Given CAMP_DIR is not a directory: ${CAMP_DIR}")
    endif()

    find_package(camp REQUIRED PATHS ${CAMP_DIR})

    message(STATUS "Checking for expected Camp target 'camp'")
    if (NOT TARGET camp)
        message(FATAL_ERROR "Camp failed to load: ${CAMP_DIR}")
    else()
        message(STATUS "Camp loaded: ${CAMP_DIR}")
        set(CAMP_FOUND TRUE CACHE BOOL "")
    endif()

    # Note: camp sets a compile feature that is not available on XL
    set_target_properties(camp PROPERTIES INTERFACE_COMPILE_FEATURES "")
else()
    message(STATUS "Camp support is OFF")
    set(CAMP_FOUND FALSE CACHE BOOL "")
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

    # Suppress warnings from cub and cuda related to the (low) version
    # of clang that XL compiler pretends to be.
    if(C_COMPILER_FAMILY_IS_XL)
        if(TARGET RAJA::cub)
            blt_add_target_definitions(
                TO RAJA::cub
                SCOPE INTERFACE
                TARGET_DEFINITIONS CUB_IGNORE_DEPRECATED_CPP_DIALECT)
        endif()

        if(TARGET cuda)
            blt_add_target_definitions(
                TO cuda
                SCOPE INTERFACE
                TARGET_DEFINITIONS THRUST_IGNORE_DEPRECATED_CPP_DIALECT)
        endif()
    endif()

else()
    message(STATUS "RAJA support is OFF" )
    set(RAJA_FOUND FALSE CACHE BOOL "")
endif()

#------------------------------------------------------------------------------
# Conduit
#------------------------------------------------------------------------------
# Find Conduit first, then find HDF5 to fix "Could NOT find HDF5" issue with
# newer CMake versions
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
# HDF5
#------------------------------------------------------------------------------
if (HDF5_DIR)
    include(cmake/thirdparty/SetupHDF5.cmake)
    blt_list_append(TO TPL_DEPS ELEMENTS hdf5)
else()
    message(STATUS "HDF5 support is OFF")
endif()

#------------------------------------------------------------------------------
# MFEM
#------------------------------------------------------------------------------
if (TARGET mfem)
    # Case: Axom included in project that also creates an mfem target, no need to recreate mfem
    # Note - white238: I can't seem to get this to pass install testing due to mfem being included
    # in multiple export sets
    message(STATUS "MFEM support is ON, using existing mfem target")
    # Add it to this export set but don't prefix it with axom::
    install(TARGETS              mfem
            EXPORT               axom-targets
            DESTINATION          lib)
    set(MFEM_FOUND TRUE CACHE BOOL "" FORCE)
elseif (MFEM_DIR)
    include(cmake/thirdparty/FindMFEM.cmake)
    # If the CMake build system was used, a CMake target for mfem already exists
    if (NOT TARGET mfem)
        # Mark mfem (and subsequent dependencies without a CMake config file) as
        # EXPORTABLE so they can be exported into axom-targets, allowing for a
        # "shrinkwrapped" CMake config
        blt_import_library( NAME       mfem
                            INCLUDES   ${MFEM_INCLUDE_DIRS}
                            LIBRARIES  ${MFEM_LIBRARIES}
                            TREAT_INCLUDES_AS_SYSTEM ON
                            EXPORTABLE ON)
        blt_list_append(TO TPL_DEPS ELEMENTS mfem)
    endif()
    blt_convert_to_system_includes(TARGET mfem)
else()
    message(STATUS "MFEM support is OFF")
endif()

# caliper-enabled mfem in device configs have extra dependencies which are not properly exported
if(TARGET mfem)
    # check if mfem depends on caliper
    set(_mfem_depends_on_caliper FALSE)
    get_target_property(_mfem_libs mfem INTERFACE_LINK_LIBRARIES)
    if("${_mfem_libs}" MATCHES "caliper")
        set(_mfem_depends_on_caliper TRUE)
    endif()

    # patch with CUDAToolkit's cupti in cuda configs
    if(_mfem_depends_on_caliper AND ENABLE_CUDA)
        if(NOT TARGET CUDA::cupti)
            find_package(CUDAToolkit REQUIRED)
        endif()

        if(TARGET CUDA::cupti)
            blt_patch_target(NAME mfem DEPENDS_ON CUDA::cupti)
        endif()
    endif()

    # patch with roctracer in hip configs
    if(_mfem_depends_on_caliper AND ENABLE_HIP)
        if(NOT TARGET roctracer)
            include(cmake/thirdparty/FindROCTracer.cmake)
            blt_list_append(TO TPL_DEPS ELEMENTS roctracer)
        endif()

        blt_patch_target(NAME mfem DEPENDS_ON roctracer)
    endif()

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
    blt_import_library( NAME       scr
                        INCLUDES   ${SCR_INCLUDE_DIRS}
                        LIBRARIES  ${SCR_LIBRARIES}
                        TREAT_INCLUDES_AS_SYSTEM ON
                        EXPORTABLE ON)
    blt_list_append(TO TPL_DEPS ELEMENTS scr)
else()
    message(STATUS "SCR support is OFF")
endif()


#------------------------------------------------------------------------------
# LUA
#------------------------------------------------------------------------------
if (LUA_DIR)
    include(cmake/thirdparty/FindLUA.cmake)
    if(NOT TARGET lua)
        blt_import_library( NAME        lua
                            LIBRARIES   ${LUA_LIBRARIES}
                            INCLUDES    ${LUA_INCLUDE_DIR}
                            EXPORTABLE  ON)
    endif()

    blt_convert_to_system_includes(TARGET lua)
    blt_list_append(TO TPL_DEPS ELEMENTS lua)
else()
    message(STATUS "LUA support is OFF")
    set(LUA_FOUND OFF CACHE BOOL "")
endif()


#------------------------------------------------------------------------------
# C2C
#------------------------------------------------------------------------------
if (C2C_DIR)
    include(cmake/thirdparty/FindC2C.cmake)
    blt_import_library(
        NAME          c2c
        INCLUDES      ${C2C_INCLUDE_DIR}
        LIBRARIES     ${C2C_LIBRARY}
        TREAT_INCLUDES_AS_SYSTEM ON
        EXPORTABLE    ON)
    blt_list_append(TO TPL_DEPS ELEMENTS c2c)
else()
    message(STATUS "c2c support is OFF")
    set(C2C_FOUND OFF CACHE BOOL "")
endif()


#------------------------------------------------------------------------------
# Remove exported OpenMP flags because they are not language agnostic
#------------------------------------------------------------------------------
set(_props)
if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13.0" )
    list(APPEND _props INTERFACE_LINK_OPTIONS)
endif()
list(APPEND _props INTERFACE_COMPILE_OPTIONS)

foreach(_target RAJA camp umpire umpire_alloc conduit::conduit)
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

# Newer versions of RAJA keeps its flags in a specific target
if(TARGET RAJA)
    get_target_property(_flags RAJA INTERFACE_LINK_LIBRARIES)
    if ( _flags )
        list(REMOVE_ITEM _flags "RAJA::openmp")
        set_target_properties( RAJA PROPERTIES INTERFACE_LINK_LIBRARIES "${_flags}" )
    endif()
endif()


#------------------------------------------------------------------------------
# Targets that need to be exported but don't have a CMake config file
#------------------------------------------------------------------------------
blt_list_append(TO TPL_DEPS ELEMENTS cuda cuda_runtime IF AXOM_ENABLE_CUDA)
blt_list_append(TO TPL_DEPS ELEMENTS blt_hip blt_hip_runtime IF AXOM_ENABLE_HIP)
blt_list_append(TO TPL_DEPS ELEMENTS openmp IF AXOM_ENABLE_OPENMP)
blt_list_append(TO TPL_DEPS ELEMENTS mpi IF AXOM_ENABLE_MPI)

foreach(dep ${TPL_DEPS})
    # If the target is EXPORTABLE, add it to the export set
    get_target_property(_is_imported ${dep} IMPORTED)
    if(NOT ${_is_imported})
        install(TARGETS              ${dep}
                EXPORT               axom-targets
                DESTINATION          lib)
        # Namespace target to avoid conflicts
        set_target_properties(${dep} PROPERTIES EXPORT_NAME axom::${dep})
    endif()
endforeach()
