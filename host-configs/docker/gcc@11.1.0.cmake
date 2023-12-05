#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/home/axom/axom_tpls/gcc-11.1.0/umpire-2023.06.0-pdygmf4s7ggj4zqbrzy436omasf6eaee;/home/axom/axom_tpls/gcc-11.1.0/raja-2023.06.0-cyz2uqefpnhonaipzcnh4oiidyrfmnsx;/home/axom/axom_tpls/gcc-11.1.0/camp-2023.06.0-csawyq2dbktfyy7srozbr7lefkz2bjmc;/home/axom/axom_tpls/gcc-11.1.0/mfem-4.5.2-dhj56qfyq3zkqiyjapqitz3iaeyjiu3r;/home/axom/axom_tpls/gcc-11.1.0/hypre-2.24.0-xsefh67wzc7yga4wotgcn4nktzsr3dko;/home/axom/axom_tpls/gcc-11.1.0/lua-5.4.4-3iwzuwuvnrymxte2nik7uqwfksk7ui5p;/home/axom/axom_tpls/gcc-11.1.0/readline-8.2-3ghsag74tthrw7kycdjwbstrhxwhiics;/home/axom/axom_tpls/gcc-11.1.0/ncurses-6.4-bqojhx5e5yk7qwezjqzchqtufyegrff7;/home/axom/axom_tpls/gcc-11.1.0/conduit-0.8.8-yohx25l7etfwzetmtcu6erdh3qnbglrn;/home/axom/axom_tpls/gcc-11.1.0/parmetis-4.0.3-wqjpma4axfmudec4ncgxkhgeyquzzneh;/home/axom/axom_tpls/gcc-11.1.0/metis-5.1.0-a6gendzgkorlnu2x6ogz4tywd4yp4c7t;/home/axom/axom_tpls/gcc-11.1.0/hdf5-1.8.22-tdsth6jhlgwhvfscnkpb24lbvvdtgobe;/home/axom/axom_tpls/gcc-11.1.0/zlib-1.2.13-k4fzfbz45b6qs3ysphaxrdowxyvjrg5p;/home/axom/axom_tpls/gcc-11.1.0/gmake-4.4.1-xoqy3bpcm7lirebrxjt4qvp3p72t3dmu;/home/axom/axom_tpls/gcc-11.1.0/blt-0.5.3-zthngkscrkhniwcjsfur5rglm6bkmge3" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@=11.1.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "-pthread" CACHE STRING "")

set(CMAKE_CXX_FLAGS "-pthread" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/bin/mpic++" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/home/axom/axom_tpls/gcc-11.1.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-yohx25l7etfwzetmtcu6erdh3qnbglrn" CACHE PATH "")

# C2C not built

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-dhj56qfyq3zkqiyjapqitz3iaeyjiu3r" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-tdsth6jhlgwhvfscnkpb24lbvvdtgobe" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-3iwzuwuvnrymxte2nik7uqwfksk7ui5p" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-cyz2uqefpnhonaipzcnh4oiidyrfmnsx" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-pdygmf4s7ggj4zqbrzy436omasf6eaee" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-csawyq2dbktfyy7srozbr7lefkz2bjmc" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


