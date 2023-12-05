#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/home/axom/axom_tpls/clang-10.0.0/umpire-2023.06.0-ffv2unko7ol7xkbswptn6reo5vkk4z36;/home/axom/axom_tpls/clang-10.0.0/raja-2023.06.0-7mjlmrtwawof7ahtizt2kh3m7pzpaa6r;/home/axom/axom_tpls/clang-10.0.0/camp-2023.06.0-cpcbuctxxtzpikmmjzt5u43oonrsqyag;/home/axom/axom_tpls/clang-10.0.0/mfem-4.5.2-vlkan3esoqo4ext3qsa6m3qiau3gcetk;/home/axom/axom_tpls/clang-10.0.0/hypre-2.24.0-a7pavarlyr7ao3qlmwxbayorpwhmoazm;/home/axom/axom_tpls/clang-10.0.0/lua-5.4.4-g77yiae3ob3q5p7bdgurcj76fgc4uxe6;/home/axom/axom_tpls/clang-10.0.0/readline-8.2-ssf54hxnbztmn4c3gekrqfmu5ovcsexe;/home/axom/axom_tpls/clang-10.0.0/ncurses-6.4-ug6qch6hz5ysjky6ebx76pzg4qumoijl;/home/axom/axom_tpls/clang-10.0.0/conduit-0.8.8-zzhc6p6w4eed6gfuu6nc6kf4uazzsy3z;/home/axom/axom_tpls/clang-10.0.0/parmetis-4.0.3-qb5ova7xrfnnjcnwjg277camap6xfu7d;/home/axom/axom_tpls/clang-10.0.0/metis-5.1.0-xlsz5xtc5npvrzlfkwjiho6uuotdmhv6;/home/axom/axom_tpls/clang-10.0.0/hdf5-1.8.22-bwa77phd2gdupwroibbaz7oa24i5s4mg;/home/axom/axom_tpls/clang-10.0.0/zlib-1.2.13-rgebxh4wloej77x6kjavjfh5cnjegtlj;/home/axom/axom_tpls/clang-10.0.0/gmake-4.4.1-ovbffagtyrtthibtg7dm4xghwjwejn2f;/home/axom/axom_tpls/clang-10.0.0/blt-0.5.3-6jop7ew6gg3noxponfcbnxswwckmjpyc" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=10.0.0
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/home/axom/axom_tpls/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "-pthread" CACHE STRING "")

set(CMAKE_CXX_FLAGS "-pthread" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,/usr/lib -Wl,-rpath,/usr/lib64" CACHE STRING "Adds a missing libstdc++ rpath")

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

set(TPL_ROOT "/home/axom/axom_tpls/clang-10.0.0" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-zzhc6p6w4eed6gfuu6nc6kf4uazzsy3z" CACHE PATH "")

# C2C not built

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-vlkan3esoqo4ext3qsa6m3qiau3gcetk" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-bwa77phd2gdupwroibbaz7oa24i5s4mg" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-g77yiae3ob3q5p7bdgurcj76fgc4uxe6" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-7mjlmrtwawof7ahtizt2kh3m7pzpaa6r" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-ffv2unko7ol7xkbswptn6reo5vkk4z36" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-cpcbuctxxtzpikmmjzt5u43oonrsqyag" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# ClangFormat disabled due to llvm and devtools not in spec

set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")

set(ENABLE_DOCS OFF CACHE BOOL "")


