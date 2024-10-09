message(STATUS "Building dependencies for Axom")
message(STATUS "CURRENT_INSTALLED_DIR -- ${CURRENT_INSTALLED_DIR}")
message(STATUS "PORT -- ${PORT}")

set(_copyright [=[
Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
other Axom Project Developers. See the top-level LICENSE file for details.

SPDX-License-Identifier: (BSD-3-Clause)
]=])

# Define a template for emitted TPL-enabled Axom host-config file
set(_host-config_hdr [=[
#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Host-config generated by vcpkg
#
# Port: @PORT@
# Architecture: @VCPKG_TARGET_ARCHITECTURE@
# Platform toolset: @VCPKG_PLATFORM_TOOLSET@
#
# vcpkg root path: @VCPKG_ROOT_DIR@
# vcpkg target triplet: @TARGET_TRIPLET@
# vcpkg target triplet file: @TARGET_TRIPLET_FILE@
#
# CMake executable path: @CMAKE_COMMAND@
#------------------------------------------------------------------------------
# To configure the code using the vcpkg toolchain:
#   cmake -C @_hc_file@  \
#         -G <generator> \
#         ../src
#
# Supported MSVC generators:
#   For x86 MSVC builds, use "Visual Studio 2022", "Visual Studio 2019" or "Visual Studio 2017"
#   For x64 MSVC builds, use "Visual Studio 2022 -A x64", "Visual Studio 2019 -A x64" or "Visual Studio 2017 Win64"
#   (note: msvc 2019 and later use the -A flag to set the architecture)
#
#
# One can also use Axom's `config-build` script:
#   cd <axom>
#   config-build.py -hc @_hc_file@                   \
#                   --msvc {2017,201764,2019,201964,2022,202264} \
#                   [other options]
#
#------------------------------------------------------------------------------
# To build the code through the command line:
#   cmake --build . --target ALL_BUILD --config Debug  [ -- -m:8 [-v:m] ]  
#
# To run tests, run either:
#   cmake --build . --target RUN_TESTS --config Debug
#   ctest -C Debug [-j8]
#
# For release builds, use 'Release' in the configuration instead of 'Debug'
#------------------------------------------------------------------------------

# Toolchain file
set(CMAKE_TOOLCHAIN_FILE "@VCPKG_ROOT_DIR@/scripts/buildsystems/vcpkg.cmake" CACHE FILEPATH "")
set(VCPKG_TARGET_TRIPLET "@TARGET_TRIPLET@" CACHE STRING "")

# CMake options
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "")

# Axom options
set(AXOM_ENABLE_TESTS ON CACHE BOOL "")
set(AXOM_ENABLE_DOCS OFF CACHE BOOL "")
set(AXOM_ENABLE_EXAMPLES ON CACHE BOOL "")

if(VCPKG_TARGET_TRIPLET MATCHES "^x64")
  set(AXOM_USE_64BIT_INDEXTYPE ON CACHE BOOL "")
endif()

# BLT options
set(ENABLE_FORTRAN OFF CACHE BOOL "")
set(ENABLE_FOLDERS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Static vs. Dynamic builds
#------------------------------------------------------------------------------
# Note: Static builds require some care and effort to get right with MSVC.  
# With a static build, choose one of
#    - disable Google Test and MSVC static MD to MT (see BLT options
#      section) or
#    - disable one of HDF5, conduit (which uses HDF5), or sidre (which
#      uses conduit).
#------------------------------------------------------------------------------

# On Windows, build shared libraries by default.
set(BUILD_SHARED_LIBS ON CACHE BOOL "")
# Shared libraries on Windows don't export symbols by default.  We'll export
# all symbols to make behavior more like Linux or Mac OS.
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON CACHE BOOL "")

# Toggle the following to disable gtest if you are compiling with static
# libraries and need HDF5
set(ENABLE_GTEST ON CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

# Toggle the following to disable changing MSVC's /MD to /MT if you are
# compiling with static libraries and need HDF5
set(BLT_ENABLE_MSVC_STATIC_MD_TO_MT ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI options
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# Set TPLs
#------------------------------------------------------------------------------
]=])

set(_conduit_dep_on [=[

set(CONDUIT_DIR "@CURRENT_INSTALLED_DIR@/share/conduit" CACHE PATH "")
set(HDF5_DIR "@CURRENT_INSTALLED_DIR@" CACHE PATH "")
]=])
set(_conduit_dep_off [=[

# conduit is disabled
# sidre requires conduit; inlet and klee require sidre
set(AXOM_ENABLE_SIDRE OFF CACHE BOOL "")
set(AXOM_ENABLE_INLET OFF CACHE BOOL "")
set(AXOM_ENABLE_KLEE OFF CACHE BOOL "")
]=])

set(_lua_dep [=[

set(LUA_DIR "@CURRENT_INSTALLED_DIR@" CACHE PATH "")
]=])

set(_mfem_dep [=[

set(MFEM_DIR "@CURRENT_INSTALLED_DIR@" CACHE PATH "")
]=])

set(_camp_dep [=[

set(CAMP_DIR "@CURRENT_INSTALLED_DIR@" CACHE PATH "")
]=])

set(_raja_dep [=[

set(RAJA_DIR "@CURRENT_INSTALLED_DIR@" CACHE PATH "")
]=])

set(_umpire_dep [=[

set(UMPIRE_DIR "@CURRENT_INSTALLED_DIR@" CACHE PATH "")
]=])


set(_opencascade_dep [=[

set(OPENCASCADE_DIR "@CURRENT_INSTALLED_DIR@" CACHE PATH "")
]=])

set(_openmp_dep [=[

# Setup OpenMP; fix MSVC linker error about unknown flag
set(ENABLE_OPENMP ON CACHE BOOL "")
set(BLT_OPENMP_LINK_FLAGS " " CACHE STRING "")
]=])

# TODO:
#  * Add features/TPLs: mpi
#  * Add tools: uncrustify, sphinx, doxygen

# Create a copyright file
file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/share/${PORT} )
set(_copyright_file ${CURRENT_PACKAGES_DIR}/share/${PORT}/copyright)
file(WRITE ${_copyright_file} "${_copyright}")

# Create a host-config file
file(MAKE_DIRECTORY ${CURRENT_PACKAGES_DIR}/include/${PORT} )
set(_hc_file ${CURRENT_PACKAGES_DIR}/include/${PORT}/hc.cmake)

# Add enabled features to host-config
message(STATUS "FEATURES: ${FEATURES}")

file(WRITE ${_hc_file}.in "${_host-config_hdr}")

if(conduit IN_LIST FEATURES)
  file(APPEND ${_hc_file}.in "${_conduit_dep_on}")
else()
  file(APPEND ${_hc_file}.in "${_conduit_dep_off}")
endif()

foreach(_dep lua mfem openmp raja umpire opencascade)
  if(${_dep} IN_LIST FEATURES)
    file(APPEND ${_hc_file}.in "${_${_dep}_dep}")
  else()
    file(APPEND ${_hc_file}.in "# ${_dep} dependency disabled")
  endif()
endforeach()

# camp is required if umpire or raja are present
if(raja IN_LIST FEATURES OR umpire IN_LIST FEATURES)
  file(APPEND ${_hc_file}.in "${_camp_dep}")
endif()

configure_file(${_hc_file}.in ${_hc_file} @ONLY)
