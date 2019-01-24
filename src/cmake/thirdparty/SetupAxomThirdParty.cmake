#------------------------------------------------------------------------------
# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

####################################
# 3rd Party Dependencies
####################################

################################
# UMPIRE
################################
if (UMPIRE_DIR)
  include(cmake/thirdparty/FindUmpire.cmake)
  blt_register_library( NAME umpire
                        INCLUDES ${UMPIRE_INCLUDE_DIRS}
                        TREAT_INCLUDES_AS_SYSTEM ON
                        LIBRARIES umpire umpire_tpl_simpool umpire_tpl_judy 
                                  umpire_resource umpire_alloc umpire_op 
                                  umpire_util umpire_strategy )
else()
  message(STATUS "Umpire support is OFF")
endif()

################################
# RAJA
################################
if (RAJA_DIR)
  include(cmake/thirdparty/FindRAJA.cmake)
  blt_register_library( NAME raja
                        INCLUDES ${RAJA_INCLUDE_DIR}
                        TREAT_INCLUDES_AS_SYSTEM ON
                        LIBRARIES ${RAJA_LIB_DIR}/libRAJA.a )
else()
  message(STATUS "RAJA support is OFF" )
endif()

################################
# Conduit
################################
if (CONDUIT_DIR)
    include(cmake/thirdparty/FindConduit.cmake)
    blt_register_library( NAME conduit
                          INCLUDES ${CONDUIT_INCLUDE_DIRS} 
                          LIBRARIES conduit
                          TREAT_INCLUDES_AS_SYSTEM ON)
    blt_register_library( NAME conduit_relay
                          INCLUDES ${CONDUIT_INCLUDE_DIRS}
                          LIBRARIES conduit_relay
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "Conduit support is OFF")
endif()


################################
# HDF5
################################
if (HDF5_DIR)
    include(cmake/thirdparty/SetupHDF5.cmake)
    blt_register_library( NAME hdf5
                          INCLUDES ${HDF5_INCLUDE_DIRS}
                          LIBRARIES ${HDF5_LIBRARIES}
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "HDF5 support is OFF")
endif()


################################
# MFEM
################################
if (MFEM_DIR)
    include(cmake/thirdparty/FindMFEM.cmake)
    blt_register_library( NAME mfem
                          INCLUDES ${MFEM_INCLUDE_DIRS}
                          LIBRARIES ${MFEM_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "MFEM support is OFF")
endif()


################################
# Shroud - Generates C/Fortran/Python bindings
################################
if(EXISTS ${SHROUD_EXECUTABLE})
    if(NOT EXISTS ${PYTHON_EXECUTABLE})
        message(FATAL_ERROR "Shroud requires PYTHON_EXECUTABLE and SHROUD_EXECUTABLE to be defined and exist.")
    endif()
    execute_process(COMMAND ${SHROUD_EXECUTABLE}
                    --cmake ${CMAKE_CURRENT_BINARY_DIR}/SetupShroud.cmake
                    ERROR_VARIABLE SHROUD_cmake_error
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(${SHROUD_cmake_error})
       message(FATAL_ERROR "Error from Shroud: ${SHROUD_cmake_error}")
    endif()
    include(${CMAKE_CURRENT_BINARY_DIR}/SetupShroud.cmake)
else()
    message(STATUS "Shroud support is OFF")
endif()


################################
# SCR
################################
if (SCR_DIR)
    include(cmake/thirdparty/FindSCR.cmake)
    blt_register_library( NAME scr
                          INCLUDES ${SCR_INCLUDE_DIRS}
                          LIBRARIES ${SCR_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON)
else()
    message(STATUS "SCR support is OFF")
endif()

