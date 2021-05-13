# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Setup HDF5
#------------------------------------------------------------------------------

# first Check for HDF5_DIR
if(NOT HDF5_DIR)
    MESSAGE(FATAL_ERROR "HDF5 support needs explicit HDF5_DIR")
endif()

# find the absolute path w/ symlinks resolved of the passed HDF5_DIR, 
# since sanity checks later need to compare against the real path
get_filename_component(HDF5_DIR_REAL "${HDF5_DIR}" REALPATH)
message(STATUS "Looking for HDF5 at: " ${HDF5_DIR_REAL})

# CMake's FindHDF5 module uses the HDF5_ROOT env var
set(HDF5_ROOT ${HDF5_DIR_REAL} CACHE PATH "" FORCE)

if(NOT WIN32)
    # use HDF5_ROOT env var for FindHDF5 with older versions of cmake
    if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
        set(ENV{HDF5_ROOT} ${HDF5_ROOT}/bin)
    endif()
endif()

# Use CMake's FindHDF5 module to locate hdf5 and setup hdf5
find_package(HDF5 REQUIRED)

# FindHDF5/find_package sets HDF5_DIR to it's installed CMake info if it exists
# we want to keep HDF5_DIR as the root dir of the install to be 
# consistent with other packages

# find the absolute path w/ symlinks resolved of the passed HDF5_DIR, 
# since sanity checks later need to compare against the real path
get_filename_component(HDF5_DIR_REAL "${HDF5_ROOT}" REALPATH)

set(HDF5_DIR ${HDF5_DIR_REAL} CACHE PATH "" FORCE)
message(STATUS "HDF5_DIR_REAL=${HDF5_DIR_REAL}")
#
# Sanity check to alert us if some how we found an hdf5 instance
# in an unexpected location.  
#
message(STATUS "Checking that found HDF5_INCLUDE_DIRS are in HDF5_DIR")

#
# HDF5_INCLUDE_DIRS may also include paths to external lib headers 
# (such as szip), so we check that *at least one* of the includes
# listed in HDF5_INCLUDE_DIRS exists in the HDF5_DIR specified.
#

# HDF5_INCLUDE_DIR is deprecated, but there are still some cases
# where HDF5_INCLUDE_DIR is set, but HDF5_INCLUDE_DIRS is not
if(NOT HDF5_INCLUDE_DIRS)
    if(HDF5_INCLUDE_DIR)
        set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "FindHDF5 did not provide HDF5_INCLUDE_DIRS or HDF5_INCLUDE_DIR.")
    endif()
endif()

if(NOT HDF5_LIBRARIES)
    message(FATAL_ERROR "FindHDF5 did not provide HDF5_LIBRARIES.")
endif()

message(STATUS "HDF5_INCLUDE_DIRS=${HDF5_INCLUDE_DIRS}")
set(check_hdf5_inc_dir_ok 0)
foreach(IDIR ${HDF5_INCLUDE_DIRS})

    # get real path of the include dir
    # w/ abs and symlinks resolved
    get_filename_component(IDIR_REAL "${IDIR}" REALPATH)
    # check if idir_real is a substring of hdf5_dir

    if("${IDIR_REAL}" MATCHES "${HDF5_DIR}")
        message(STATUS " ${IDIR_REAL} includes HDF5_DIR (${HDF5_DIR})")
        set(check_hdf5_inc_dir_ok 1)
    endif()
endforeach()

if(NOT check_hdf5_inc_dir_ok)
    message(FATAL_ERROR " ${HDF5_INCLUDE_DIRS} does not include HDF5_DIR")
endif()

#
# filter HDF5_LIBRARIES to remove hdf5_hl if it exists
# we don't use hdf5_hl, but if we link with it will become
# a transitive dependency
#
set(HDF5_HL_LIB FALSE)
foreach(LIB ${HDF5_LIBRARIES})
    if("${LIB}" MATCHES "hdf5_hl")
        set(HDF5_HL_LIB ${LIB})
    endif()
endforeach()

if(HDF5_HL_LIB)
    message(STATUS "Removing hdf5_hl from HDF5_LIBRARIES")
    list(REMOVE_ITEM HDF5_LIBRARIES ${HDF5_HL_LIB})
endif()


#
# Display main hdf5 cmake vars
#
message(STATUS "HDF5 Include Dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 Libraries:    ${HDF5_LIBRARIES}")
message(STATUS "HDF5 Definitions:  ${HDF5_DEFINITIONS}")
message(STATUS "HDF5 is parallel:  ${HDF5_IS_PARALLEL}")

# if newer style hdf5 imported targets exist, use those on windows
if(WIN32 AND TARGET hdf5::hdf5-shared AND BUILD_SHARED_LIBS)
    # reg shared ver of imported lib target
    message(STATUS "HDF5 using hdf5::hdf5-shared target")
    blt_import_library(NAME      hdf5
                       LIBRARIES hdf5::hdf5-shared
                       INCLUDES  ${HDF5_INCLUDE_DIRS}
                       TREAT_INCLUDES_AS_SYSTEM ON
                       EXPORTABLE ON)
elseif(WIN32 AND TARGET hdf5::hdf5-static )
    # reg static ver of imported lib target
    message(STATUS "HDF5 using hdf5::hdf5-static target")
    blt_import_library(NAME      hdf5
                       LIBRARIES hdf5::hdf5-static
                       EXPORTABLE ON)
elseif(TARGET hdf5)
    # legacy hdf5 CMake build system support creates an hdf5 target we use directly
    message(STATUS "HDF5 using hdf5 target")

    set_property(TARGET hdf5 
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 "${HDF5_INCLUDE_DIRS}")
else()
    # reg includes and libs with blt
    message(STATUS "HDF5 using HDF5_DEFINITIONS + HDF5_INCLUDE_DIRS + HDF5_LIBRARIES")
    message(STATUS "HDF5_DEFINITIONS:  ${HDF5_DEFINITIONS}")
    message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5_LIBRARIES:    ${HDF5_LIBRARIES}")
    blt_import_library(NAME hdf5
                       DEFINES   ${HDF5_DEFINITIONS}
                       INCLUDES  ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARIES}
                       TREAT_INCLUDES_AS_SYSTEM ON
                       EXPORTABLE ON )
endif()
