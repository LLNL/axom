#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=16.0.2
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/rocmcc/rocmcc-5.6.0beta1-magic/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/rocmcc/rocmcc-5.6.0beta1-magic/bin/amdclang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/rocmcc/rocmcc-5.6.0beta1-magic/bin/amdflang" CACHE PATH "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_05_15_44_11/clang-16.0.2" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-dkkcbcxe2pzcqr7rywpdyufu4grrouko" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-t4lvw3ilqlupgl2wsl66solflqacb227" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0-vzzp4owphc7c44j7ilrpp3vf2vqpcfjq" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-irckv3y42x74zonmwzzkw7f6cgzksf7a" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-jrfceuovnaha3kyv4qfnilsvvmorpvfa" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-33zoxgfc2z4czau6tq3qepdi4es3nhst" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-ruzaqttp57mpklysyn7nrb7znx6re5b2" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-m2fdamzxdaucg7xc4wrs3qtzibfc33av" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


