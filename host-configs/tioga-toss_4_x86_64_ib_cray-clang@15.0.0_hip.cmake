#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=15.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/opt/rocm-5.4.3/llvm/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/opt/rocm-5.4.3/llvm/bin/amdclang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/opt/rocm-5.4.3/llvm/bin/amdflang" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-Mfreeform" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3/bin/mpif90" CACHE PATH "")

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

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags -L/opt/rocm-5.4.3/hip/../llvm/lib -L/opt/rocm-5.4.3/hip/lib -Wl,-rpath,/opt/rocm-5.4.3/hip/../llvm/lib:/opt/rocm-5.4.3/hip/lib -lpgmath -lflang -lflangrti -lompstub -lamdhip64  -L/opt/rocm-5.4.3/hip/../lib64 -Wl,-rpath,/opt/rocm-5.4.3/hip/../lib64 " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_06_30_13_06_22/clang-15.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-khnnbowe63malwwrqpi2jmq45os276ra" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-pxcld4oajhkqj4rpdfpm4hk6ep3avwwv" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0-5eugp2rnn6gsmqyxtg4o6oupdinkqbig" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-t3m2n7h4nosb6yo2hf2rkwl2bsutqeed" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-7rzfp3hkjm43x2hwssloyenublslxg6h" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-v7wpbzo6h5uctubr5yrrgwxkfyhvgdqp" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-vq2jxjgp2xzqgdshwxtit2ngxl5wfozi" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-wlqmpth747wftntclnpzofqpdv7jyhjk" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


