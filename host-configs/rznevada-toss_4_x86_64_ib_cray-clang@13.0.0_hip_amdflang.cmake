#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@13.0.0_hip_amdflang
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray_temp_dir/toss_4_x86_64_ib_cray/2022_01_28_08_04_25/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray_temp_dir/toss_4_x86_64_ib_cray/2022_01_28_08_04_25/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray_temp_dir/toss_4_x86_64_ib_cray/2022_01_28_08_04_25/spack/lib/spack/env/clang/flang" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/rocmcc-tce/rocmcc-4.3.1/bin/amdclang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/rocmcc-tce/rocmcc-4.3.1/bin/hipcc" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/rocmcc-tce/rocmcc-4.3.1/bin/amdflang" CACHE PATH "")

endif()

set(CMAKE_Fortran_FLAGS "-Mfreeform" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.9-rocmcc-4.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.9-rocmcc-4.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.9-rocmcc-4.3.1/bin/mpif90" CACHE PATH "")

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

set(HIP_ROOT_DIR "/opt/rocm-4.3.1/hip" CACHE STRING "")

set(CMAKE_HIP_ARCHITECTURES "gfx908" CACHE STRING "")

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--disable-new-dtags -L/opt/rocm-4.3.1/hip/../llvm/lib -L/opt/rocm-4.3.1/hip/lib -Wl,-rpath,/opt/rocm-4.3.1/hip/../llvm/lib:/opt/rocm-4.3.1/hip/lib -lpgmath -lflang -lflangrti -lompstub -lamdhip64 -L/opt/rocm-4.3.1/hip/../lib64 -Wl,-rpath,/opt/rocm-4.3.1/hip/../lib64 -lhsakmt -lamd_comgr" CACHE STRING "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

# Root directory for generated TPLs

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib_cray_temp_dir/toss_4_x86_64_ib_cray/2022_01_28_08_04_25/clang-13.0.0_hip_amdflang" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.7.2axom" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.3.0" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.2.0" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.3.5" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-0.14.0" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-6.0.0" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_3_x86_64_ib/2020_08_21_22_18_57/gcc-8.1.0" CACHE PATH "")

# ClangFormat disabled due to disabled devtools

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/python3.7" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-1.87/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


