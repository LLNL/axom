#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/umpire-2023.06.0-bftoz5wfuvfwldjwe3szsi54hdpw2hfi;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/raja-2023.06.0-tv34emklwgpl4ddrggbpkthu6utweo2a;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/camp-2023.06.0-vuk2xmgpzxl2xb2jupm3uxqxbgrzrmbf;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/lua-5.4.4-e7cdhkbaxuyqljoyozzxrt2n57erhszx;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/conduit-0.8.8-3li2jcc2mecqbanoew5tznr4uzl2k2jg;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/parmetis-4.0.3-c2ub4mddqt7faykccls5eae7xogpi32w;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/metis-5.1.0-t6qk2khpceidsormndy6gj56qbnehxai;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/hdf5-1.8.22-nydokvym7klek6hqoyijbwf3tbmd32ux;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/zlib-1.2.13-zd2x25bmmjwxhw7siucfywjgk7oc262j;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/c2c-1.8.0-zgwmbota3vhdxcdjlwxbhomlwllwhvnc;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1/blt-0.5.3-zghdcshuqlqnbrezcylwnefie7pkvsab;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/usr/tce/packages/clang/clang-10.0.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/graphviz-7.1.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/doxygen-1.8.14;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/cppcheck-2.9;/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1;/usr/tcetmp;/usr/tce/packages/cmake/cmake-3.21.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@=8.3.1.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE STRING "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE "/usr/tce/packages/gcc/gcc-4.9.3/lib64;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3" CACHE STRING "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_26_18_11_07/gcc-8.3.1.1" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-3li2jcc2mecqbanoew5tznr4uzl2k2jg" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-zgwmbota3vhdxcdjlwxbhomlwllwhvnc" CACHE PATH "")

# MFEM not built

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-nydokvym7klek6hqoyijbwf3tbmd32ux" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-e7cdhkbaxuyqljoyozzxrt2n57erhszx" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-tv34emklwgpl4ddrggbpkthu6utweo2a" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-bftoz5wfuvfwldjwe3szsi54hdpw2hfi" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-vuk2xmgpzxl2xb2jupm3uxqxbgrzrmbf" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/2023_04_18_13_41_48/._view/srxt35kojgk77f2222mk6mgv7z5jyyzz" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


