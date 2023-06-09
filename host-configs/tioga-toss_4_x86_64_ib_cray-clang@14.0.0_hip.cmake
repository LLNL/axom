#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=14.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_05_19_09_36_48/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_05_19_09_36_48/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_05_19_09_36_48/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/opt/rocm-5.2.3/llvm/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/opt/rocm-5.2.3/llvm/bin/amdclang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/opt/rocm-5.2.3/llvm/bin/amdflang" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-Mfreeform" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.16-rocmcc-5.2.3/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.16-rocmcc-5.2.3/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.16-rocmcc-5.2.3/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# HIP

#------------------------------------------------------------------------------


set(ENABLE_HIP ON CACHE BOOL "")

set(HIP_ROOT_DIR "/opt/rocm-5.2.3/hip" CACHE STRING "")

set(HIP_CLANG_INCLUDE_PATH "/opt/rocm-5.2.3/hip/../llvm/lib/clang/14.0.0/include" CACHE PATH "")

set(CMAKE_CXX_FLAGS "--std=c++14" CACHE STRING "")

set(CMAKE_HIP_ARCHITECTURES "gfx90a" CACHE STRING "")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags -L/opt/rocm-5.2.3/hip/../llvm/lib -L/opt/rocm-5.2.3/hip/lib -Wl,-rpath,/opt/rocm-5.2.3/hip/../llvm/lib:/opt/rocm-5.2.3/hip/lib -lpgmath -lflang -lflangrti -lompstub -lamdhip64  -L/opt/rocm-5.2.3/hip/../lib64 -Wl,-rpath,/opt/rocm-5.2.3/hip/../lib64 -lhsakmt " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_05_19_09_36_48/clang-14.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-7e7dezqy3yd4xniycht5jsusfunmjghg" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-jhjgx666nyzstsoypqlb3nltxfldr2yh" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.0-7wiw2z3ffyoeegiqaz46b5bqmlox42sd" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-63q4hszz6x5m5bpme25u4bmwjgucxkt5" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-oboalbf3obehbhpxwiijtadx7gra2l6t" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-5magqvbudhw35i7hnwi2qbfs73lf76pw" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-i2nfecyspzw7hmmmzvdlkdwadakrh2ri" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-vatpalnxoikyyieg37wuspflxkauldzh" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


