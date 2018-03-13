##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-clang@4.0.0
##################################

# cmake from uberenv
# cmake executable path: /collab/usr/global/tools/cmake/bgqos_0/cmake-3.8.2/bin/cmake

#######
# using clang@4.0.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set(ENABLE_FORTRAN OFF CACHE BOOL "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_11_13_15_16_50/spack/opt/spack/bgqos_0/clang-4.0.0" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-y7l7rg3y7j7hxdsom3nwcpwwfq3tu75b" CACHE PATH "")

# scr not built by uberenv

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-lps3i452m3y4ebespcl7atv4bukml7xq" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3.2-qtbnjv44avhlwb7l6zl3taypjzllki4w" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-22q34tajrlunpnme5ezp72t4iwy2qpgz" CACHE PATH "")

# python not built by uberenv

# lua not built by uberenv

# doxygen not built by uberenv

# sphinx not built by uberenv

# shroud not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc bgq clang@4.0.0 host configs
##############################################################################

set(ENABLE_DOCS    OFF CACHE BOOL "")
set(ENABLE_PYTHON  OFF CACHE BOOL "")

set(CMAKE_SKIP_RPATH TRUE CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI      ON CACHE BOOL "")
# Note: On BGQ, CMake uses wrong linker flags when using find MPI, This allows us use the wrapper directly via
#   the CMake's compiler variables.
set(ENABLE_FIND_MPI OFF CACHE BOOL "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC TRUE CACHE BOOL "Ensures that tests will be wrapped with srun to run on the backend nodes")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

