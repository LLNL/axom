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
# host-config for 'sqa-uno' Windows machine-- using Visual Studio 15
#
# Run the following from a build dir for a 32-bit configuration
#   cmake -G "Visual Studio 15 2017"                       \
#         -C ..\host-configs\other\sqa-uno-MSVC.cmake      \
#         <path-to-axom>
#
# Run the following from a build dir for a 64-bit configuration
#   cmake -G "Visual Studio 15 2017 Win64"               \
#         -C ..\host-configs\other\sqa-uno-MSVC.cmake    \
#         <path-to-axom> 
#
# On Tangelo: 
#   /c/Users/yeh14/Program/cmake-3.11.2-win64-x64/bin/cmake.exe -G "Visual Studio 15 2017 Win64" -C ../../Axom/host-configs/other/tangelo.cmake ../../Axom/
#
# Build the code from the command line as follows (/m for parallel build in msbuild):
#   cmake --build . --config {Release, Debug, RelWithDebInfo} [-- //m:8]
#
# Test the code as follows (j for parallel testing):
#   ctest -j8 -C {Release,Debug,RelWithDebInfo}
# 
# Install the come from the command line as follows:
#   cmake --build . --config Release --target install
#
#------------------------------------------------------------------------------

# Set the HOME variable (%USERPROFILE% in Windows)
string(REPLACE "\\" "/" HOME "$ENV{USERPROFILE}")

### Setup some devtool/TPL paths

# Enable sidre using conduit but not hdf5 (until we resolve configuration problems with hdf5)
set(ENABLE_SIDRE OFF CACHE BOOL "")
#set(HDF5_DIR    ... CACHE PATH "")
#set(CONDUIT_DIR "${HOME}/Projects/install/conduit-no-hdf5" CACHE PATH "")

# Note: Doxygen assumes graphviz 'dot' is in PATH
#set(DOXYGEN_EXECUTABLE "${HOME}/Chocolatey/bin/doxygen.exe" CACHE PATH "") 

# Note: uncrustify disabled since this version is different than uberenv version and changes whitespace
# set(UNCRUSTIFY_EXECUTABLE "${HOME}/Code/UniversalIndentGUI/indenters/uncrustify.exe" CACHE PATH "")

# Setup MPI
set(ENABLE_MPI OFF CACHE BOOL "")
#set(MPI_HOME  "C:/Program Files/Microsoft HPC Pack 2008 R2" CACHE PATH "")
#set(MPI_GUESS_LIBRARY_NAME "MSMPI" CACHE STRING "")


### Set some additional options
set(ENABLE_FOLDERS ON CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#making the error go away
set(BLT_CXX_FLAGS "-D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING=1" CACHE STRING "")
