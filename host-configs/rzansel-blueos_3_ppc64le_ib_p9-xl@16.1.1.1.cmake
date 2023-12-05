#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/umpire-2023.06.0-cmose4tqzlwsrjjhoklstt6ga3yym6ni;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/raja-2023.06.0-uqr6nybua6ic5xjwqd75hqtngiz6wjwb;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/camp-2023.06.0-3k5jzo2qageywkgo7i65ybqksjqbu67j;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/mfem-4.5.2-5tpknd5eztbbomidbr3qcs75hyvdqegw;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/hypre-2.24.0-kgmrxqyc3pjlwgzvgbchceub5u4v4tgt;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/lua-5.4.4-22nibaxxfzwid5jk3udivaptba2xtw4r;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/conduit-0.8.8-tnyxf6f7sh7kg3pzvt5acfp4btsprsp7;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/parmetis-4.0.3-zsfvlzuxwixta5mjt4dqwikslhx6l6dd;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/metis-5.1.0-hpfbnba3mgk2opmqn3efr2mmda76yacg;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/hdf5-1.8.22-eduvgrxspycy6xfw2t6i6uogdtzgnygm;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/zlib-1.2.13-tnpzx5nxltbrlx27lozjehkg2lkb7lew;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/c2c-1.8.0-vxyvicuyswybrekljpj7khd56uobpc62;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1/blt-0.5.3-zleps37vlw57krmv2zruxq6bpxpn6fik;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/usr/tcetmp/packages/lapack/lapack-3.9.0-P9-xl-2020.11.12;/usr/tce/packages/clang/clang-10.0.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/graphviz-7.1.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/doxygen-1.9.6;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/cppcheck-2.9;/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19;/usr/tcetmp;/usr/tce/packages/cmake/cmake-3.21.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: xl@=16.1.1.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/spack/lib/spack/env/xl/xlc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/spack/lib/spack/env/xl/xlc++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/spack/lib/spack/env/xl/xlf90" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/xl/xl-2022.08.19/bin/xlc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/xl/xl-2022.08.19/bin/xlC" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/xl/xl-2022.08.19/bin/xlf2003" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_CXX_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_Fortran_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19/bin/mpixlc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19/bin/mpixlC" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19/bin/mpixlf" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE STRING "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS "${BLT_EXE_LINKER_FLAGS} -Wl,-rpath,/usr/tce/packages/xl/xl-2022.08.19/lib" CACHE STRING "Adds a missing rpath for libraries associated with the fortran compiler")

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath,/usr/tce/packages/xl/xl-2022.08.19/lib" CACHE STRING "Adds a missing rpath for libraries associated with the fortran compiler")

set(BLT_FORTRAN_FLAGS "-WF,-C!  -qxlf2003=polymorphic" CACHE STRING "Converts C-style comments to Fortran style in preprocessed files")

set(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE "/usr/tce/packages/gcc/gcc-4.9.3/lib64;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3" CACHE STRING "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-tnyxf6f7sh7kg3pzvt5acfp4btsprsp7" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-vxyvicuyswybrekljpj7khd56uobpc62" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-5tpknd5eztbbomidbr3qcs75hyvdqegw" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-eduvgrxspycy6xfw2t6i6uogdtzgnygm" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-22nibaxxfzwid5jk3udivaptba2xtw4r" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-uqr6nybua6ic5xjwqd75hqtngiz6wjwb" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-cmose4tqzlwsrjjhoklstt6ga3yym6ni" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-3k5jzo2qageywkgo7i65ybqksjqbu67j" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/2023_10_17_10_41_04/._view/5gug4et2qymmckk7yvtub4qcapp3figm" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


