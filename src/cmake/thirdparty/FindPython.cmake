###############################################################################
# Copyright (c) 2014-2015, Lawrence Livermore National Security, LLC.
# 
# Produced at the Lawrence Livermore National Laboratory
# 
# LLNL-CODE-666778
# 
# All rights reserved.
# 
# This file is part of Conduit. 
# 
# For details, see: http://llnl.github.io/conduit/.
# 
# Please also read conduit/LICENSE
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the disclaimer below.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
# 
###############################################################################

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
        COMMAND PYTHONPATH=${CMAKE_Python_MODULE_DIRECTORY}
            ${PYTHON_EXECUTABLE} setup.py -v
            build
            install
              --install-purelib=${CMAKE_Python_MODULE_DIRECTORY}
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

##------------------------------------------------------------------------------
## blt_add_python_module
##
## Creates a shared library to be used as a Python module.
## All options to blt_add_library may be used.
##
## The library is created in CMAKE_Python_MODULE_DIRECTORY
##------------------------------------------------------------------------------
macro(blt_add_python_module)
    include_directories(${PYTHON_INCLUDE_DIR})

    # Force shard libraries
    blt_add_library(${ARGV} SHARED)

    # Python wants the name to be 'name.so', without leading 'lib'
    set_target_properties(${arg_NAME}
        PROPERTIES
          PREFIX ""
          LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Python_MODULE_DIRECTORY}
    )
endmacro(blt_add_python_module)
