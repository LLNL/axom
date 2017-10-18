##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# bgqos_0-clang@3.7.0
##################################

# cmake from uberenv
# cmake executable path: /usr/local/tools/cmake-3.4.3/bin/cmake

#######
# using clang@3.7.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/collab/usr/gapps/opnsrc/gnu/dev/lnx-2.12-ppc/bgclang/r284961-stable/llnl/bin/mpiclang++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set(ENABLE_FORTRAN OFF CACHE BOOL "")

# Root directory for generated TPLs
set(TPL_ROOT "/usr/workspace/wsrzc/axom/thirdparty_libs/builds/2017_07_19_14_04_28/spack/opt/spack/bgqos_0/clang-3.7.0" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.16-bosaqxj3xd5fhyovqnda3rgj2kjsj4ah" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.2.1-oiiieme5mlcpao7pqwrk2mdquxnaguqm" CACHE PATH "")

# mfem from uberenv
set(MFEM_DIR "${TPL_ROOT}/mfem-3.3-2ctejypbka5bdvs43cnv3twskptqpfjs" CACHE PATH "")

# boost headers from uberenv
set(BOOST_DIR "${TPL_ROOT}/boost-headers-1.58.0-qddl3bajxtossmhy4mazvjpah4zgx5aj" CACHE PATH "")

# python not built by uberenv

# lua not built by uberenv

# doxygen not built by uberenv

# sphinx not built by uberenv

# uncrustify not built by uberenv

# lcov and genhtml not built by uberenv

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc bgq clang@3.7.0 host configs
##############################################################################

set(ENABLE_DOCS    OFF CACHE BOOL "")
set(ENABLE_PYTHON  OFF CACHE BOOL "")

set(CMAKE_SKIP_RPATH TRUE CACHE BOOL "")

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE BOOL "")
set(ENABLE_FIND_MPI OFF CACHE BOOL "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC TRUE CACHE BOOL "Ensures that tests will be wrapped with srun to run on the backend nodes")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

