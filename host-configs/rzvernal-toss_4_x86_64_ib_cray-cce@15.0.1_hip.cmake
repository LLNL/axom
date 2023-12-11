#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/umpire-2023.06.0-znpux7jjwwb3gb2hr4vi6nbhyrkfo5mo;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/raja-2023.06.0-ovzqujoxykloaaitzbtkg6egqp6qjr2v;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/camp-2023.06.0-nsepimnez7uuv3dxqiji2twih3swa6cr;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/mfem-4.5.2-7rrla6up5iou2fic4riu776n4isvt4wf;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/hypre-2.24.0-ci2qpkdifjovd7ki22dscgl6goyjefai;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/lua-5.4.4-gduannh2dcxphs7odrtkrpnmabgcikom;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/ncurses-6.4-cnq7lljngizvsdryq3csevc3tbh4k7vh;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/conduit-0.8.8-w5ybmdn2qvfbkooexazuzxnmvyc4zfyl;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/parmetis-4.0.3-u2udk5n4552zcftdccg5ol33mtulfq2i;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/metis-5.1.0-lssxpp454fp6c4fbcl2d6rwugzcrpgvm;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/hdf5-1.8.22-uvl6nsaahk6vl3mvoq6a432lontoo6zv;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/c2c-1.8.0-drh2lcdrfn627tcfkocetim7snaakicw;/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1/blt-0.5.3-xxkf276amm4zj5wm5dx44clspivj42vx;/opt/rocm-5.4.3;/opt/rocm-5.4.3/llvm;/opt/rocm-5.4.3;/opt/rocm-5.4.3/hip;/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.25-rocmcc-5.4.3-cce-15.0.1;/usr/tce" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: cce@=15.0.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/spack/lib/spack/env/cce/craycc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/spack/lib/spack/env/cce/case-insensitive/crayCC" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/spack/lib/spack/env/cce/crayftn" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/cce-tce/cce-15.0.1/bin/craycc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/cce-tce/cce-15.0.1/bin/crayCC" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/cce-tce/cce-15.0.1/bin/crayftn" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-ef" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -DNDEBUG" CACHE STRING "")

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

set(BLT_CMAKE_IMPLICIT_LINK_LIBRARIES_EXCLUDE "unwind" CACHE STRING "")

set(CMAKE_EXE_LINKER_FLAGS " -L/opt/rocm-5.4.3/hip/../lib64 -Wl,-rpath,/opt/rocm-5.4.3/hip/../lib64  -L/opt/rocm-5.4.3/hip/../lib -Wl,-rpath,/opt/rocm-5.4.3/hip/../lib -lamd_comgr -lhsa-runtime64 " CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray/2023_07_26_18_12_45/cce-15.0.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-w5ybmdn2qvfbkooexazuzxnmvyc4zfyl" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-drh2lcdrfn627tcfkocetim7snaakicw" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-7rrla6up5iou2fic4riu776n4isvt4wf" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-uvl6nsaahk6vl3mvoq6a432lontoo6zv" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-gduannh2dcxphs7odrtkrpnmabgcikom" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-ovzqujoxykloaaitzbtkg6egqp6qjr2v" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-znpux7jjwwb3gb2hr4vi6nbhyrkfo5mo" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-nsepimnez7uuv3dxqiji2twih3swa6cr" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


