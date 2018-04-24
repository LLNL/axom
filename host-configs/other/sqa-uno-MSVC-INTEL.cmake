#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------
# host-config for 'sqa-uno' Windows machine
# using the Intel 18 toolchain for Visual Studio 15
#
# Run the following from a build dir for a 64-bit configuration:
#   cmake -G "Visual Studio 15 2017 Win64"                   \
#         -T "Intel C++ Compiler 18.0"                       \
#         -C ..\host-configs\other\sqa-uno-MSVC-intel.cmake  \
#         <path-to-axom>
#
# Build the code from the command line as follows (/m for parallel build in msbuild):
#   cmake --build . --config {Release, Debug, RelWithDebInfo} [-- /m:8]
#
# Test the code as follows (j for parallel testing):
#   ctest -j8 -C {Release,Debug,RelWithDebInfo}
# 
# Install the come from the command line as follows:
#   cmake --build . --config Release --target install
#
# Note: MPI in this configuration requires an initial login.
# If the MPI tests hang, try to run the blt_mpi_smoke test manually and enter 
# your credentials. They should work automatically (e.g. through ctest) after that.
#------------------------------------------------------------------------------

### Setup Axom components

# Disable Sidre until we can successfully link to HDF5 and Conduit on Windows
set(ENABLE_SIDRE OFF CACHE BOOL "")
#set(HDF5_DIR    ... CACHE PATH "")
#set(CONDUIT_DIR ... CACHE PATH "")

# Disable Fortran since Visual Studio is not generating modules for our component Fortan interfaces
set(ENABLE_FORTRAN OFF CACHE BOOL "")


### Setup some devtool/TPL paths

# Set the HOME variable (%USERPROFILE% in Windows)
string(REPLACE "\\" "/" HOME "$ENV{USERPROFILE}")

set(BOOST_DIR "${HOME}/Code/boost" CACHE PATH "")

# Note: Doxygen assumes graphviz 'dot' is in PATH
set(DOXYGEN_EXECUTABLE "${HOME}/Chocolatey/bin/doxygen.exe" CACHE PATH "") 

# Note: uncrustify disabled since this version is different than uberenv version and changes whitespace
# set(UNCRUSTIFY_EXECUTABLE "${HOME}/Code/UniversalIndentGUI/indenters/uncrustify.exe" CACHE PATH "")

# Setup MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_HOME              "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2018.1.156/windows/mpi/intel64" CACHE PATH "")
set(MPI_C_COMPILER        "${MPI_HOME}/bin/mpiicc.bat" CACHE PATH "")
set(MPI_CXX_COMPILER      "${MPI_HOME}/bin/mpiicpc.bat" CACHE PATH "")
set(MPI_Fortran_COMPILER  "${MPI_HOME}/bin/mpiifort.bat" CACHE PATH "")


### Set some additional options
set(ENABLE_FOLDERS ON CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")
