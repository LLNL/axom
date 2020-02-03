# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# 3rd Party Dependencies
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# UMPIRE
#------------------------------------------------------------------------------
if (UMPIRE_DIR)
    include(cmake/thirdparty/FindUmpire.cmake)
    blt_register_library( NAME      umpire
                          INCLUDES  ${UMPIRE_INCLUDE_DIRS}
                          LIBRARIES umpire
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "Umpire support is OFF")
endif()


#------------------------------------------------------------------------------
# RAJA
#------------------------------------------------------------------------------
if (RAJA_DIR)
    include(cmake/thirdparty/FindRAJA.cmake)
else()
    message(STATUS "RAJA support is OFF" )
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
# HDF5
#------------------------------------------------------------------------------
if (HDF5_DIR)
    include(cmake/thirdparty/SetupHDF5.cmake)
    blt_register_library( NAME      hdf5
                          INCLUDES  ${HDF5_INCLUDE_DIRS}
                          LIBRARIES ${HDF5_LIBRARIES}
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "HDF5 support is OFF")
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
    if(NOT EXISTS ${PYTHON_EXECUTABLE})
        message(FATAL_ERROR "Shroud requires PYTHON_EXECUTABLE and SHROUD_EXECUTABLE to be defined and exist.")
    endif()
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
foreach(_target RAJA camp umpire umpire_alloc)
    if(TARGET ${_target})
        get_target_property(_flags ${_target} INTERFACE_COMPILE_OPTIONS)
        if ( _flags )

            message(STATUS "Removing OpenMP Flags from target[${_target}]")

            string( REPLACE "${OpenMP_CXX_FLAGS}" ""
                    correct_flags "${_flags}" )
            string( REPLACE "${OpenMP_Fortran_FLAGS}" ""
                    correct_flags "${correct_flags}" )

            set_target_properties( ${_target} PROPERTIES INTERFACE_COMPILE_OPTIONS
                                   "${correct_flags}" )
        endif()
    endif()
endforeach()
