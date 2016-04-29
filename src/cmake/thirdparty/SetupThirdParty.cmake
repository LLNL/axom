####################################
# 3rd Party Dependencies
####################################

################################
# Conduit
################################
if (CONDUIT_DIR)
  include(cmake/thirdparty/FindConduit.cmake)
  blt_register_library( NAME conduit
                        INCLUDES ${CONDUIT_INCLUDE_DIRS} )
  blt_register_library( NAME conduit_io
                        INCLUDES ${CONDUIT_INCLUDE_DIRS} )
endif()


################################
# HDF5
################################
if (HDF5_DIR)
  include(cmake/thirdparty/FindHDF5.cmake)
  blt_register_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARY} )
endif()


################################
# Sparsehash
################################
if (SPARSEHASH_DIR)
  include(cmake/thirdparty/FindSparsehash.cmake)
  blt_register_library(NAME sparsehash
                       INCLUDES ${SPARSEHASH_INCLUDE_DIRS})
endif()


################################
# Documentation Packages
################################
if (DOXYGEN_EXECUTABLE)
  find_package(Doxygen)
endif()

if (SPHINX_EXECUTABLE)
  include(cmake/thirdparty/FindSphinx.cmake)
endif()


################################
# linting via Uncrustify
################################
if (UNCRUSTIFY_EXECUTABLE)
  include(cmake/thirdparty/FindUncrustify.cmake)
endif()


################################
# Find boost headers
################################
if (ENABLE_BOOST)
  if (DEFINED BOOST_ROOT)
    find_package(Boost
                 1.55
                 REQUIRED)
    blt_register_library(NAME boost
                         INCLUDES ${Boost_INCLUDE_DIR})
    MESSAGE(STATUS "Boost include dir: " ${Boost_INCLUDE_DIR})
    MESSAGE(STATUS "Boost version: " ${Boost_VERSION})
  else()
    MESSAGE(FATAL_ERROR "ENABLE_BOOST is true, but BOOST_ROOT was not set.  Check your host-config file.")
  endif()
endif()

################################
# Python
################################

if(ENABLE_PYTHON AND PYTHON_EXECUTABLE)
    ################################
    # Setup includes for Python & Numpy
    ################################
    include(FindPython)
    message(STATUS "Using Python Include: ${PYTHON_INCLUDE_DIRS}")
    include_directories(${PYTHON_INCLUDE_DIRS})
    # if we don't find python, throw a fatal error
    if(NOT PYTHON_FOUND)
        message(FATAL_ERROR "ENABLE_EXECUTABLE is set, but Python wasn't found.")
    endif()

    ## Set the Python module directory
    # relative path (used with install)
    set(CMAKE_Python_SITE_PACKAGES
        "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages"
        CACHE PATH
        "Relative path where all Python modules will go in the install tree"
    )
    # build site-packages
    set(CMAKE_Python_MODULE_DIRECTORY
        "${PROJECT_BINARY_DIR}/${CMAKE_Python_SITE_PACKAGES}"
        CACHE PATH
        "Directory where all Python modules will go in the build tree"
    )

    file(MAKE_DIRECTORY ${CMAKE_Python_MODULE_DIRECTORY})
    set(ENV{PYTHONPATH} ${CMAKE_Python_MODULE_DIRECTORY})

    INSTALL(DIRECTORY DESTINATION ${CMAKE_Python_SITE_PACKAGES})
    INSTALL(CODE " set(ENV\{PYTHONPATH\} ${CMAKE_INSTALL_PREFIX}/${CMAKE_Python_SITE_PACKAGES}) ")
endif(ENABLE_PYTHON AND PYTHON_EXECUTABLE)

################################
# Lua
################################

if (LUA_DIR)
    # Set the hint for FindLua
    set (ENV{LUA_DIR}  ${LUA_DIR})
    find_package(Lua)

    set(CMAKE_Lua_MODULE_DIRECTORY
        "${PROJECT_BINARY_DIR}/lib/lua"
        CACHE PATH
        "Directory where all Lua modules will go in the build tree"
    )

    set(LUA_EXECUTABLE ${LUA_DIR}/bin/lua)
endif()

##------------------------------------------------------------------------------
## blt_add_lua_module
##
## Creates a shared library to be used as a Lua module.
## All options to blt_add_library may be used.
##
## The library is created in CMAKE_Lua_MODULE_DIRECTORY
## by default. OUTPUT_DIR can be used to change the location.
##
## NAME is the name of the Lua module.
## The target name will be ${arg_NAME}-lua-module and the
## library is named ${arg_NAME}.so.
## This allow lib${arg_NAME}.a to also be created by using
## blt_add_library directly.
## 
##------------------------------------------------------------------------------
macro(blt_add_lua_module)
    include_directories(${LUA_INCLUDE_DIR})

    set(singleValueArgs NAME )
    ## parse the arguments
    ## only parse NAME, blt_add_library will do the real work
    cmake_parse_arguments(arg_module
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

    # Force shard libraries
    blt_add_library(
        OUTPUT_DIR ${CMAKE_Lua_MODULE_DIRECTORY}
        ${ARGV}
        NAME ${arg_module_NAME}-lua-module
        SHARED
    )

    # Lua wants the name to be 'name.so', without leading 'lib'
    set_target_properties(${arg_module_NAME}-lua-module PROPERTIES
        PREFIX ""
	OUTPUT_NAME ${arg_module_NAME}
    )
endmacro(blt_add_lua_module)
