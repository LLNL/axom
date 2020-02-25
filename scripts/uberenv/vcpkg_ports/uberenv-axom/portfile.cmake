if(VCPKG_CMAKE_SYSTEM_NAME STREQUAL "WindowsStore")
    message(FATAL_ERROR "${PORT} does not currently support UWP")
endif()

include(vcpkg_common_functions)

message(STATUS "Building dependencies for Axom")

message(STATUS "CURRENT_INSTALLED_DIR -- ${CURRENT_INSTALLED_DIR}")
message(STATUS "PORT -- ${PORT}")

set(_copyright [=[
Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
other Axom Project Developers. See the top-level COPYRIGHT file for details.

SPDX-License-Identifier: (BSD-3-Clause)
]=])

# Define a template for emitted TPL-enabled Axom host-config file
set(_host-config_hdr [=[
#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Host-config generated by vcpkg
#
# Port: @PORT@
# Architecture: @VCPKG_TARGET_ARCHITECTURE@
# Platform toolset: @VCPKG_PLATFORM_TOOLSET@
#
# vcpkg root path: @VCPKG_ROOT_PATH@
# vcpkg target triplet: @TARGET_TRIPLET@
# vcpkg target triplet file: @TARGET_TRIPLET_FILE@
#
# CMake system name: @CMAKE_SYSTEM_NAME@
# CMake system version: @CMAKE_SYSTEM_VERSION@
# CMake executable path: @CMAKE_COMMAND@
#------------------------------------------------------------------------------
# Empty/useless variables:
#   VS path: @VCPKG_VISUAL_STUDIO_PATH@
#   VC Package root: @VCPKG_ROOT@
#   Linkage: @VCPKG_CRT_LINKAGE@
#   Library linkage: @VCPKG_CRT_LINKAGE@
#   CMake system name: @VCPKG_CMAKE_SYSTEM_NAME@
#   CMake system version: @VCPKG_CMAKE_SYSTEM_VERSION@
#------------------------------------------------------------------------------
# To configure the code using the vcpkg toolchain:
#   cmake -C @_hc_file@ ../src
#
# To build the code through the command line:
#   cmake --build . --target ALL_BUILD --config Debug  [ -- -m:8 [-v:m] ]  
#
# To run tests, run either:
#   cmake --build . --target RUN_TESTS --config Debug
#   ctest -C Debug [-j8]
#------------------------------------------------------------------------------

# On Windows, build shared libraries by default.
set(BUILD_SHARED_LIBS ON CACHE BOOL "")
# Static builds require some care and effort to get right.  With a static
# build, choose one of
#    - disable Google Test and MSVC static MD to MT (see BLT options
#      section) or
#    - disable one of HDF5, conduit (which uses HDF5), or sidre (which
#      uses conduit).

# Toolchain file
set(CMAKE_TOOLCHAIN_FILE @VCPKG_ROOT_PATH@/scripts/buildsystems/vcpkg.cmake CACHE FILEPATH "")
set(VCPKG_TARGET_TRIPLET @TARGET_TRIPLET@ CACHE STRING "")

# Set TPLs
set(CONDUIT_DIR @CURRENT_INSTALLED_DIR@/share/conduit CACHE PATH "")
set(HDF5_DIR @CURRENT_INSTALLED_DIR@ CACHE PATH "")

# Axom options
set(AXOM_ENABLE_TESTS ON CACHE BOOL "")
set(AXOM_ENABLE_DOCS OFF CACHE BOOL "")
set(AXOM_ENABLE_EXAMPLES ON CACHE BOOL "")
# set(AXOM_ENABLE_SIDRE OFF CACHE BOOL "")

# BLT options
set(ENABLE_FORTRAN OFF CACHE BOOL "")
set(ENABLE_FOLDERS ON CACHE BOOL "")
# Toggle the following to disable gtest if you are compiling with static
# libraries and need HDF5
set(ENABLE_GTEST ON CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")
# Toggle the following to disable changing MSVC's /MD to /MT if you are
# compiling with static libraries and need HDF5
set(BLT_ENABLE_MSVC_STATIC_MD_TO_MT ON CACHE BOOL "")

# MPI options
set(ENABLE_MPI OFF CACHE BOOL "")
# If MSMPI and no other MPI is installed, turn ENABLE_MPI ON and CMake
# should automatically detect it.  If CMake doesn't auto-detect MSMPI,
# or if you need to use another MPI, you will need to specify the MPI
# compiler wrappers and helper settings.
#
# Here are settings that might be appropriate for Intel compiler:
#
# set(MPI_C_COMPILER "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpicc.bat" CACHE PATH "")
# set(MPI_CXX_COMPILER "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpicc.bat" CACHE PATH "")
# set(MPI_Fortran_COMPILER "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpifc.bat" CACHE PATH "")
# set(MPIEXEC "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpiexec.exe" CACHE PATH "")
# set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")
#
# Here are example settings pointing to MSMPI (use when Intel is installed):
#
# set(MPI_C_COMPILER "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpicc.bat" CACHE PATH "")
# set(MPI_CXX_COMPILER "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpicc.bat" CACHE PATH "")
# set(MPI_Fortran_COMPILER "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpifc.bat" CACHE PATH "")
# set(MPIEXEC "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.5.281/windows/mpi/intel64/bin/mpiexec.exe" CACHE PATH "")
# set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

# cmake options
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "")

# TODO:
#  * Add TPLs: mfem, umpire, raja
#  * Add tools: uncrustify, sphinx, doxygen

# DONE:
#  * Add conduit with HDF5
#  * Add hints to get MPI working
#  * Add vcpkg toolchain file -- CMAKE_TOOLCHAIN_FILE
#  * Set vcpkg triplet -- VCPKG_TARGET_TRIPLET 

]=])

# Create a copyright file
file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/share/${PORT} )
set(_copyright_file ${CURRENT_PACKAGES_DIR}/share/${PORT}/copyright)
file(WRITE ${_copyright_file} "${_copyright}")

# Create a host-config file
file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/include/${PORT} )
set(_hc_file ${CURRENT_PACKAGES_DIR}/include/${PORT}/hc.cmake)

file(WRITE ${_hc_file}.in ${_host-config_hdr})
configure_file(${_hc_file}.in ${_hc_file} @ONLY)

