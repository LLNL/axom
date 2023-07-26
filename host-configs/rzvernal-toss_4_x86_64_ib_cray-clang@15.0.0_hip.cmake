#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/umpire-2022.10.0-ddf5mqkyfusgjme7tmexdlldnplwxnh4;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/raja-2022.10.5-shf5smaswum4au5aehbyuteaoun3wgsu;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/camp-2022.10.1-jf4xeumelu43v7l2oqgm3o2loigbrnqe;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/mfem-4.5.2-qu3gj4kd7xsetoefz23bofifgyufllhi;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/hypre-2.24.0-qbm2npjkwsp4vryixh65wbi2vohvpvbv;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/lua-5.4.4-nciboohbuioaqb35wwibitzc3ypfjbvp;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/ncurses-6.4-aqfvdlf2fuogw7pvfsg33qvvsdskwxn3;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/conduit-0.8.8-hqffdgfu3fhgmjlctcuywp66hp5ldg7e;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/parmetis-4.0.3-22wxjkvf6g6l74z2ztuk3imydnmpi7gq;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/metis-5.1.0-i4fcjmmgu4qquorg5pjbuipwd7m3viy6;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/hdf5-1.8.22-eoxieadxuddhem2eot4dqfdj47xngqpj;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/c2c-1.8.0-da5avt3ikppaabh4do3vgximikvvsjhy;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0/blt-0.5.3-57fuyfjl4ksfli3zskcbvykyrwtrwfne;/opt/rocm-5.4.3;/opt/rocm-5.4.3/llvm;/opt/rocm-5.4.3;/opt/rocm-5.4.3/hip;/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3;/usr/tce" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=15.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/opt/rocm-5.4.3/llvm/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/opt/rocm-5.4.3/llvm/bin/amdclang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/opt/rocm-5.4.3/llvm/bin/amdflang" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-Mfreeform" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

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

#------------------------------------------------
# ROCm
#------------------------------------------------

set(HIP_ROOT_DIR "/opt/rocm-5.4.3/hip" CACHE PATH "")

set(HIP_CXX_COMPILER "/opt/rocm-5.4.3/hip/bin/hipcc" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx90a" CACHE STRING "")

set(AMDGPU_TARGETS "gfx90a" CACHE STRING "")

set(GPU_TARGETS "gfx90a" CACHE STRING "")

#------------------------------------------------------------------------------

# Axom ROCm specifics

#------------------------------------------------------------------------------


set(ENABLE_HIP ON CACHE BOOL "")

set(HIP_CLANG_INCLUDE_PATH "/opt/rocm-5.4.3/hip/../llvm/lib/clang/15.0.0/include" CACHE PATH "")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags -L/opt/rocm-5.4.3/hip/../llvm/lib -L/opt/rocm-5.4.3/hip/lib -Wl,-rpath,/opt/rocm-5.4.3/hip/../llvm/lib:/opt/rocm-5.4.3/hip/lib -lpgmath -lflang -lflangrti -lompstub -lamdhip64  -L/opt/rocm-5.4.3/hip/../lib64 -Wl,-rpath,/opt/rocm-5.4.3/hip/../lib64  -L/opt/rocm-5.4.3/hip/../lib -Wl,-rpath,/opt/rocm-5.4.3/hip/../lib -lamd_comgr -lhsa-runtime64 " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_13_16_52_45/clang-15.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-hqffdgfu3fhgmjlctcuywp66hp5ldg7e" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-da5avt3ikppaabh4do3vgximikvvsjhy" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-qu3gj4kd7xsetoefz23bofifgyufllhi" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-eoxieadxuddhem2eot4dqfdj47xngqpj" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-nciboohbuioaqb35wwibitzc3ypfjbvp" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.10.5-shf5smaswum4au5aehbyuteaoun3wgsu" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.10.0-ddf5mqkyfusgjme7tmexdlldnplwxnh4" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-jf4xeumelu43v7l2oqgm3o2loigbrnqe" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


