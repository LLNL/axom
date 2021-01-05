#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# SYS_TYPE: linux-ubuntu16.04-ivybridge
# Compiler Spec: gcc@8.1.0
#------------------------------------------------------------------------------
# CMake executable path: /usr/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "")

set(CMAKE_C_FLAGS "-pthread" CACHE PATH "")

set(CMAKE_CXX_FLAGS "-pthread" CACHE PATH "")

set(BLT_CXX_STD "c++14" CACHE PATH "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

# Root directory for generated TPLs
set(TPL_ROOT "/home/axom/axom_tpls/gcc-8.1.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.6.0" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.1.0" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.21" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.3.5" CACHE PATH "")

# SCR not built

set(RAJA_DIR "${TPL_ROOT}/raja-0.12.1" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-4.0.1" CACHE PATH "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/bin/mpic++" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE PATH "")

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(ENABLE_DOCS OFF CACHE BOOL "")

# ClangFormat disabled due to disabled devtools
set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

#------------------------------------------------------------------------------
# Other machine specifics
#------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")


