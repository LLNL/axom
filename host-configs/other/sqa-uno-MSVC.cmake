#------------------------------------------------------------------------------
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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
# host-config for 'sqa-uno' Windows machine-- using visual studio 15
# Run the following from a build dir to configure for 32- and 64-bit configurations:
#   cmake -G "Visual Studio 15 2017"       -C ..\host-configs\other\sqa-uno-MSVC.cmake <path-to-axom>
#   cmake -G "Visual Studio 15 2017 Win64" -C ..\host-configs\other\sqa-uno-MSVC.cmake <path-to-axom>
#------------------------------------------------------------------------------

## Setup Axom components
# Disable Sidre until we can successfully link to HDF5 and Conduit on Windows
set(ENABLE_SIDRE OFF CACHE BOOL "")
#set(HDF5_DIR    ... CACHE PATH "")
#set(CONDUIT_DIR ... CACHE PATH "")

## Setup some devtool/TPL paths
# Set the HOME variable -- it is %USERPROFILE% in Windows
string(REPLACE "\\" "/" HOME "$ENV{USERPROFILE}")
set(BOOST_DIR "${HOME}/Code/boost" CACHE PATH "")
# assumes graphviz 'dot' is in PATH
set(DOXYGEN_EXECUTABLE "${HOME}/Chocolatey/bin/doxygen.exe" CACHE PATH "") 
# Note: uncrustify disabled since this version is different than uberenv version and changes whitespace
# set(UNCRUSTIFY_EXECUTABLE "${HOME}/Code/UniversalIndentGUI/indenters/uncrustify.exe" CACHE PATH "")

## Setup MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_GUESS_LIBRARY_NAME "MSMPI" CACHE STRING "")
set(MPI_HOME  "C:/Program Files/Microsoft HPC Pack 2008 R2" CACHE PATH "")

## Set some additional options
set(ENABLE_FOLDERS ON CACHE BOOL "")
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")
