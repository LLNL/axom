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

####################################
# 3rd Party Dependencies
####################################

################################
# Conduit
################################
if (CONDUIT_DIR)
    include(cmake/thirdparty/FindConduit.cmake)
    blt_register_library( NAME conduit
                          INCLUDES ${CONDUIT_INCLUDE_DIRS} 
                          LIBRARIES  conduit)
    blt_register_library( NAME conduit_relay
                          INCLUDES ${CONDUIT_INCLUDE_DIRS}
                          LIBRARIES  conduit_relay)
else()
    message(STATUS "Conduit support is OFF")
endif()

################################
# HDF5
################################
if (HDF5_DIR)
    include(cmake/thirdparty/SetupHDF5.cmake)
    blt_register_library(NAME hdf5
                         INCLUDES ${HDF5_INCLUDE_DIRS}
                         LIBRARIES ${HDF5_LIBRARIES} )
else()
    message(STATUS "HDF5 support is OFF")
endif()

################################
# MFEM
################################
if (MFEM_DIR)
    include(cmake/thirdparty/FindMFEM.cmake)
    blt_register_library( NAME mfem
                          INCLUDES  ${MFEM_INCLUDE_DIRS}
                          LIBRARIES ${MFEM_LIBRARY} )
else()
    message(STATUS "MFEM support is OFF")
endif()


################################
# Setup toolkit generate targets
################################

if(EXISTS ${SHROUD_EXECUTABLE})
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
# Python
################################

if(ENABLE_PYTHON AND PYTHON_EXECUTABLE)
    ################################
    # Setup includes for Python
    ################################
    include(cmake/thirdparty/FindPython.cmake)
    message(STATUS "Using Python Include: ${PYTHON_INCLUDE_DIRS}")
    # if we don't find python, throw a fatal error
    if(NOT PYTHON_FOUND)
        message(FATAL_ERROR "ENABLE_PYTHON is set, but Python wasn't found.")
    endif()

    ## Set the Python module directory
    # relative path (used with install)
    set(BLT_Python_SITE_PACKAGES
        "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
        CACHE PATH
        "Relative path where all Python modules will go in the install tree"
    )
    # build site-packages
    set(BLT_Python_MODULE_DIRECTORY
        "${PROJECT_BINARY_DIR}/${BLT_Python_SITE_PACKAGES}"
        CACHE PATH
        "Directory where all Python modules will go in the build tree"
    )

    file(MAKE_DIRECTORY ${BLT_Python_MODULE_DIRECTORY})
    set(ENV{PYTHONPATH} ${BLT_Python_MODULE_DIRECTORY})

    INSTALL(DIRECTORY DESTINATION ${BLT_Python_SITE_PACKAGES})
    INSTALL(CODE " set(ENV\{PYTHONPATH\} ${CMAKE_INSTALL_PREFIX}/${BLT_Python_SITE_PACKAGES}) ")

    blt_register_library(
        NAME python
        INCLUDES ${PYTHON_INCLUDE_DIRS}
        LIBRARIES ${PYTHON_LIBRARIES}
    )
else()
    message(STATUS "Python support is OFF")
endif(ENABLE_PYTHON AND PYTHON_EXECUTABLE)

################################
# Lua
################################

if (LUA_DIR)
    # Set the hint for FindLua
    set (ENV{LUA_DIR}  ${LUA_DIR})
    find_package(Lua)
    if(NOT LUA_FOUND)
        message( FATAL_ERROR "Requested Lua was not found: ${LUA_DIR}")
    endif()

    set(LUA_EXECUTABLE ${LUA_DIR}/bin/lua)

    blt_register_library(
        NAME lua
        INCLUDES ${LUA_INCLUDE_DIR}
        LIBRARIES ${LUA_LIBRARIES}
    )
else()
    message(STATUS "LUA support is OFF")
endif()

################################
# SCR
################################
if (SCR_DIR)
    include(cmake/thirdparty/FindSCR.cmake)
    blt_register_library(NAME scr
                         INCLUDES ${SCR_INCLUDE_DIRS}
                         LIBRARIES ${SCR_LIBRARY} )
else()
    message(STATUS "SCR support is OFF")
endif()


