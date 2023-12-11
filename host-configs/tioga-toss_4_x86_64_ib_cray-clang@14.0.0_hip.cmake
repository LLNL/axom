#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/umpire-2023.06.0-5shsoyik6cl37unfukgypuy5wq3xrf5d;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/raja-2023.06.0-mihlyyxgy5ksjvblaztoa7lernlfuq5c;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/camp-2023.06.0-kjnl6w2c7g7373mmyy6z2kmnvwqcjiat;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/mfem-4.5.2-z5scxsmafi5skftjzrszn2cgfg5smgxr;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/hypre-2.24.0-kirtgmu2tqj5qeewzcgwxekiw2mjjvnq;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/lua-5.4.4-s6rzmz2h2k53z53grnpe237uwwzr4nz3;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/ncurses-6.4-wwijdtdjzvmf4224vzbzcntupmetijri;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/conduit-0.8.8-2kuebxwabtkkcikcqclqytfllwm3qwhg;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/parmetis-4.0.3-g4rtt7uakthoyf32g4ug6nntsxsh3im3;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/metis-5.1.0-bhi2eebvzv5maenn2vbmxsr2ftlqqv5n;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/hdf5-1.8.22-fo5wthv72ys3e263kvxykonmixg2sxun;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/c2c-1.8.0-56n4sm5x4ecvwmdazemuqvikg5ngwdwf;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0/blt-0.5.3-2z65cg5imrjg7ktgft4ob5tgxg2cvn4f;/opt/rocm-5.2.3;/opt/rocm-5.2.3/llvm;/opt/rocm-5.2.3;/opt/rocm-5.2.3/hip;/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.16-rocmcc-5.2.3;/usr/tce" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=14.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/opt/rocm-5.2.3/llvm/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/opt/rocm-5.2.3/llvm/bin/amdclang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/opt/rocm-5.2.3/llvm/bin/amdflang" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-Mfreeform" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.16-rocmcc-5.2.3/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.16-rocmcc-5.2.3/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.16-rocmcc-5.2.3/bin/mpif90" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(MPIEXEC_EXECUTABLE "/usr/global/tools/flux_wrappers/bin/srun" CACHE PATH "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# ROCm
#------------------------------------------------

set(HIP_ROOT_DIR "/opt/rocm-5.2.3/hip" CACHE PATH "")

set(HIP_CXX_COMPILER "/opt/rocm-5.2.3/hip/bin/hipcc" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx90a" CACHE STRING "")

set(AMDGPU_TARGETS "gfx90a" CACHE STRING "")

set(GPU_TARGETS "gfx90a" CACHE STRING "")

#------------------------------------------------------------------------------

# Axom ROCm specifics

#------------------------------------------------------------------------------


set(ENABLE_HIP ON CACHE BOOL "")

set(HIP_CLANG_INCLUDE_PATH "/opt/rocm-5.2.3/hip/../llvm/lib/clang/14.0.0/include" CACHE PATH "")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags -L/opt/rocm-5.2.3/hip/../llvm/lib -L/opt/rocm-5.2.3/hip/lib -Wl,-rpath,/opt/rocm-5.2.3/hip/../llvm/lib:/opt/rocm-5.2.3/hip/lib -lpgmath -lflang -lflangrti -lompstub -lamdhip64  -L/opt/rocm-5.2.3/hip/../lib64 -Wl,-rpath,/opt/rocm-5.2.3/hip/../lib64  -L/opt/rocm-5.2.3/hip/../lib -Wl,-rpath,/opt/rocm-5.2.3/hip/../lib -lamd_comgr -lhsa-runtime64 " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_17_41_54/clang-14.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-2kuebxwabtkkcikcqclqytfllwm3qwhg" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-56n4sm5x4ecvwmdazemuqvikg5ngwdwf" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-z5scxsmafi5skftjzrszn2cgfg5smgxr" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-fo5wthv72ys3e263kvxykonmixg2sxun" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-s6rzmz2h2k53z53grnpe237uwwzr4nz3" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-mihlyyxgy5ksjvblaztoa7lernlfuq5c" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-5shsoyik6cl37unfukgypuy5wq3xrf5d" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-kjnl6w2c7g7373mmyy6z2kmnvwqcjiat" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


