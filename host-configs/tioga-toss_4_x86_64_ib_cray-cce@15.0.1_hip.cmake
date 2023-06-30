#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: cce@=15.0.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/spack/lib/spack/env/cce/craycc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/spack/lib/spack/env/cce/case-insensitive/crayCC" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/spack/lib/spack/env/cce/crayftn" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/cce-tce/cce-15.0.1/bin/craycc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/cce-tce/cce-15.0.1/bin/crayCC" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/cce-tce/cce-15.0.1/bin/crayftn" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-ef" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_CXX_FLAGS "-O1" CACHE STRING "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3-cce-15.0.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3-cce-15.0.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3-cce-15.0.1/bin/mpif90" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(MPIEXEC_EXECUTABLE "/usr/global/tools/flux_wrappers/bin/srun" CACHE PATH "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# HIP

#------------------------------------------------------------------------------


set(ENABLE_HIP ON CACHE BOOL "")

set(HIP_ROOT_DIR "/opt/rocm-5.4.3/hip" CACHE STRING "")

set(HIP_CLANG_INCLUDE_PATH "/opt/rocm-5.4.3/hip/../llvm/lib/clang/15.0.0/include" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx90a" CACHE STRING "")

set(BLT_CMAKE_IMPLICIT_LINK_LIBRARIES_EXCLUDE "unwind" CACHE STRING "")

set(CMAKE_EXE_LINKER_FLAGS " -L/opt/rocm-5.4.3/hip/../lib64 -Wl,-rpath,/opt/rocm-5.4.3/hip/../lib64 " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/cce-15.0.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-3vqsbwogrjipvqva7asa4eoth4xhsem7" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-wk54ms543nc2nssjqaqfknyrqynvglcd" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0-v6blj6vtfsp6mdzqa2vfaf6hjt4blamf" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-kyuaayh7js6xvo4ylxtj2toquq5ixzve" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-hfobxinovtz5o2a23ibcmifuwyjc52ha" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-r72afubbarpoq6dxaibcpkxmgra4awpo" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-msxhtslz6cdqkz3zmhpanfbzvqp47ju3" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-z6i3l6j5xgqoryscu4bx3b65da3ubndn" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


