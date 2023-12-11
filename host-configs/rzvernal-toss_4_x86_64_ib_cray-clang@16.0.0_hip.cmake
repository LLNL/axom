#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/umpire-2023.06.0-yeccv6q7zasapkqmulwukmay2xfq5hyu;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/raja-2023.06.0-7y7neucz73mp5aidhw7kfbmqhwgsr4ww;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/camp-2023.06.0-zu25serllpaiepsclsd3o7w4clj5k2bs;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/mfem-4.5.2-znyrzghqsshldarhqa6gweqgukahiyjw;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/hypre-2.24.0-hodcjhfvx7c6fupyaw5kgoiuazbxcaij;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/lua-5.4.4-j446jxtxyu4x2byxkcvzv4ek4w5pr3fh;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/ncurses-6.4-sqvzpcbeunczs72a55mzrqi7rpl7aojx;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/conduit-0.8.8-nexjzr5p6zy3lb7kbkojaqbmvn64bwnf;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/parmetis-4.0.3-l3yx62tarv37nnultafaw2q7yows6c22;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/metis-5.1.0-gwyy4vhrzd77gbdxd4thjxtxwrjads4l;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/hdf5-1.8.22-gph5gvans7w7rhb3acho7c3yydyxbck6;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/c2c-1.8.0-oht7wdi5u5r4zlf7mcdk36xcuvt3y5j7;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0/blt-0.5.3-no74mmcw3sf324mpvpdzlu57kmg2xwby;/opt/rocm-5.6.0;/opt/rocm-5.6.0/llvm;/opt/rocm-5.6.0;/opt/rocm-5.6.0/hip;/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.6.0;/usr/tce" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=16.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/opt/rocm-5.6.0/llvm/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/opt/rocm-5.6.0/llvm/bin/amdclang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/opt/rocm-5.6.0/llvm/bin/amdflang" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-Mfreeform" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -DNDEBUG" CACHE STRING "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.6.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.6.0/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.6.0/bin/mpif90" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(MPIEXEC_EXECUTABLE "/usr/global/tools/flux_wrappers/bin/srun" CACHE PATH "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# ROCm
#------------------------------------------------

set(HIP_ROOT_DIR "/opt/rocm-5.6.0/hip" CACHE PATH "")

set(HIP_CXX_COMPILER "/opt/rocm-5.6.0/hip/bin/hipcc" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx90a" CACHE STRING "")

set(AMDGPU_TARGETS "gfx90a" CACHE STRING "")

set(GPU_TARGETS "gfx90a" CACHE STRING "")

#------------------------------------------------------------------------------

# Axom ROCm specifics

#------------------------------------------------------------------------------


set(ENABLE_HIP ON CACHE BOOL "")

set(HIP_CLANG_INCLUDE_PATH "/opt/rocm-5.6.0/hip/../llvm/lib/clang/16.0.0/include" CACHE PATH "")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags -L/opt/rocm-5.6.0/hip/../llvm/lib -L/opt/rocm-5.6.0/hip/lib -Wl,-rpath,/opt/rocm-5.6.0/hip/../llvm/lib:/opt/rocm-5.6.0/hip/lib -lpgmath -lflang -lflangrti -lompstub -lamdhip64  -L/opt/rocm-5.6.0/hip/../lib64 -Wl,-rpath,/opt/rocm-5.6.0/hip/../lib64  -L/opt/rocm-5.6.0/hip/../lib -Wl,-rpath,/opt/rocm-5.6.0/hip/../lib -lamd_comgr -lhsa-runtime64 " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/clang-16.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-nexjzr5p6zy3lb7kbkojaqbmvn64bwnf" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-oht7wdi5u5r4zlf7mcdk36xcuvt3y5j7" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-znyrzghqsshldarhqa6gweqgukahiyjw" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-gph5gvans7w7rhb3acho7c3yydyxbck6" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-j446jxtxyu4x2byxkcvzv4ek4w5pr3fh" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-7y7neucz73mp5aidhw7kfbmqhwgsr4ww" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-yeccv6q7zasapkqmulwukmay2xfq5hyu" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-zu25serllpaiepsclsd3o7w4clj5k2bs" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


