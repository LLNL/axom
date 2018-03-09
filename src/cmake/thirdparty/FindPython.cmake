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


# Find the interpreter first
if(PYTHON_DIR AND NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE ${PYTHON_DIR}/bin/python)
endif()

find_package(PythonInterp REQUIRED)
if(PYTHONINTERP_FOUND)
        
        MESSAGE(STATUS "PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE}")
        
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
                                "import sys;from distutils.sysconfig import get_python_inc;sys.stdout.write(get_python_inc())"
                        OUTPUT_VARIABLE PYTHON_INCLUDE_DIR
                        ERROR_VARIABLE ERROR_FINDING_INCLUDES)
        MESSAGE(STATUS "PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR}")

        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
                                "import sys;from distutils.sysconfig import get_python_lib;sys.stdout.write(get_python_lib())"
                        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_DIR
                        ERROR_VARIABLE ERROR_FINDING_SITE_PACKAGES_DIR)
        MESSAGE(STATUS "PYTHON_SITE_PACKAGES_DIR ${PYTHON_SITE_PACKAGES_DIR}")

        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
                                "import sys;from distutils.sysconfig import get_config_var; sys.stdout.write(get_config_var('LIBDIR'))"
                        OUTPUT_VARIABLE PYTHON_LIB_DIR
                        ERROR_VARIABLE ERROR_FINDING_LIB_DIR)
        MESSAGE(STATUS "PYTHON_LIB_DIR ${PYTHON_LIB_DIR}")
        
        # check for python libs differs for windows python installs
        if(NOT WIN32)
            set(PYTHON_GLOB_TEST "${PYTHON_LIB_DIR}/libpython*")
        else()
            set(PYTHON_GLOB_TEST "${PYTHON_LIB_DIR}/python*.lib")
        endif()
            
        FILE(GLOB PYTHON_GLOB_RESULT ${PYTHON_GLOB_TEST})
        get_filename_component(PYTHON_LIBRARY "${PYTHON_GLOB_RESULT}" ABSOLUTE)
        MESSAGE(STATUS "{PythonLibs from PythonInterp} using: PYTHON_LIBRARY=${PYTHON_LIBRARY}")
        find_package(PythonLibs)
        
        if(NOT PYTHONLIBS_FOUND)
            MESSAGE(FATAL_ERROR "Failed to find Python Libraries using PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")
        endif()
        
endif()


find_package_handle_standard_args(Python  DEFAULT_MSG
                                  PYTHON_LIBRARY PYTHON_INCLUDE_DIR)


FUNCTION(PYTHON_ADD_DISTUTILS_SETUP target_name)
    MESSAGE(STATUS "Configuring python distutils setup: ${target_name}")
    set(timestamp ${CMAKE_CURRENT_BINARY_DIR}/${target_name}.time)
    add_custom_command(
        OUTPUT  ${timestamp}
        COMMAND PYTHONPATH=${BLT_Python_MODULE_DIRECTORY}
            ${PYTHON_EXECUTABLE} setup.py -v
            build
            install
              --install-purelib=${BLT_Python_MODULE_DIRECTORY}
              --install-scripts=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            COMMAND ${CMAKE_COMMAND} -E touch ${timestamp}
            DEPENDS  ${setup_file} ${ARGN}
    )

    add_custom_target(${target_name} ALL DEPENDS ${timestamp})

    # also use distutils for the install ...
    INSTALL(CODE "
        EXECUTE_PROCESS(
            COMMAND
                ${PYTHON_EXECUTABLE} setup.py -v
                build
                install
                   --prefix=${CMAKE_INSTALL_PREFIX}
            OUTPUT_VARIABLE PY_DIST_UTILS_INSTALL_OUT
            ERROR_VARIABLE PY_DIST_UTILS_INSTALL_ERR
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        MESSAGE(STATUS \"Installing: ${target_name}\\n\${PY_DIST_UTILS_INSTALL_OUT}\")
        if(PY_DIST_UTILS_INSTALL_ERR)
            MESSAGE(STATUS \"ERR  \${PY_DIST_UTILS_INSTALL_ERR}\")
        endif()
        ")

ENDFUNCTION(PYTHON_ADD_DISTUTILS_SETUP)

#FUNCTION(PYTHON_ADD_HYBRID_MODULE target_name dest_dir py_name setup_file py_sources)
#    MESSAGE(STATUS "Configuring hybrid python module: ${target_name}")
#    PYTHON_ADD_DISTUTILS_SETUP("${target_name}_py_setup"
#                               ${dest_dir}
#                               ${setup_file}
#                               ${py_sources})
#    PYTHON_ADD_MODULE(${target_name} ${ARGN})
#    SET_TARGET_PROPERTIES(${target_name} PROPERTIES
#                                         LIBRARY_OUTPUT_DIRECTORY   
#                             ${CMAKE_BINARY_DIR}/${dest_dir}/${py_name})
#
#    install(TARGETS ${target_name}
#            EXPORT  conduit
#            LIBRARY DESTINATION ${dest_dir}/${py_name}
#            ARCHIVE DESTINATION ${dest_dir}/${py_name}
#            RUNTIME DESTINATION ${dest_dir}/${py_name}
#    )
#
#ENDFUNCTION(PYTHON_ADD_HYBRID_MODULE)
