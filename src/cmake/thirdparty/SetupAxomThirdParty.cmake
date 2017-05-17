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
endif()

################################
# HDF5
################################
if (HDF5_DIR)
  include(cmake/thirdparty/SetupHDF5.cmake)
  blt_register_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARIES} )
endif()

################################
# MFEM
################################
if (MFEM_DIR)
  include(cmake/thirdparty/FindMFEM.cmake)
  blt_register_library( NAME mfem
                        INCLUDES  ${MFEM_INCLUDE_DIRS}
                        LIBRARIES ${MFEM_LIBRARY} )
endif()

################################
# Find boost headers
################################
if (BOOST_DIR)
    include(cmake/thirdparty/SetupBoost.cmake)
    blt_register_library(NAME boost
                         INCLUDES ${Boost_INCLUDE_DIR})
endif()

################################
# Setup toolkit generate targets
################################

if(EXISTS ${SHROUD_EXECUTABLE})
    execute_process(COMMAND ${SHROUD_EXECUTABLE} --sitedir
                    OUTPUT_VARIABLE SHROUD_Site
                    ERROR_VARIABLE SHROUD_Site_Error
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(${SHROUD_Site_Error})
       message(FATAL_ERROR "Error from Shroud: ${SHROUD_Site_Error}")
    endif()
    include(${SHROUD_Site}/cmake/SetupShroud.cmake)
endif()

macro(uncrustify_shroud)
    # Must be after a call to add_shroud macro
    # Only run uncrustify if shroud has just run
    # XXX uncrustify.cfg is hardwired
    set(_cfg ${PROJECT_SOURCE_DIR}/uncrustify.cfg)
    set(_uncrustify ${CMAKE_CURRENT_BINARY_DIR}/${_basename}.uncrustify)
    if(UNCRUSTIFY_FOUND AND (EXISTS ${_cfg}))
        add_custom_command(
            OUTPUT ${_uncrustify}
            DEPENDS  ${_timestamp}
            COMMAND ${UNCRUSTIFY_EXECUTABLE}
                    -c ${_cfg} --no-backup `cat ${_cfiles}`
            COMMAND touch ${_uncrustify}
            COMMENT "Running uncrustify for ${arg_YAML_INPUT_FILE}."
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        add_custom_target(${_shroud_target}_uncrustify
            DEPENDS ${_uncrustify}
        )
        add_dependencies(generate ${_shroud_target}_uncrustify)
    endif()
endmacro(uncrustify_shroud)

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
endif()

