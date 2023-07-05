#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=16.0.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/rocm/rocm-5.6.0beta1/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/rocm/rocm-5.6.0beta1/bin/amdclang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/rocm/rocm-5.6.0beta1/bin/amdflang" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-Mfreeform" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_CXX_FLAGS "-O1" CACHE STRING "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich/cray-mpich-8.1.25-rocmcc-5.6.0beta1-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich/cray-mpich-8.1.25-rocmcc-5.6.0beta1-magic/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich/cray-mpich-8.1.25-rocmcc-5.6.0beta1-magic/bin/mpif90" CACHE PATH "")

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

set(HIP_ROOT_DIR "/usr/tce/packages/rocm/rocm-5.6.0beta1/hip" CACHE STRING "")

set(HIP_CLANG_INCLUDE_PATH "/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/../llvm/lib/clang/16.0.0/include" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx90a" CACHE STRING "")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags -L/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/../llvm/lib -L/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/lib -Wl,-rpath,/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/../llvm/lib:/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/lib -lpgmath -lflang -lflangrti -lompstub -lamdhip64  -L/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/../lib64 -Wl,-rpath,/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/../lib64  -L/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/../lib -Wl,-rpath,/usr/tce/packages/rocm/rocm-5.6.0beta1/hip/../lib " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/clang-16.0.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-7adubjhv5u42sm3d4twybv6qt7gsvqpn" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-kr75tuh775f4c3pdw4bvspfcxj3uyqhz" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0-dxjexcmgoslsmpilwbrnetdwb72zyy6a" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-asogc3ke4zd4wb45az6seiukwt2u2gak" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-5wtwr4wd7ypvtgnpgx3t5rh5dto7n4hr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-ss52llt7w6p7o6oc62io5aw6ufgpkuvl" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-y376pnonynkm2frlpliibht3zxm3lxru" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-pqnbwjhc3aezfebrgldia26dr2gbnudj" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


