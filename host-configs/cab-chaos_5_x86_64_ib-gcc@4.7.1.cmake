##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-gcc@4.7.1
##################################

# cmake from uberenv
# cmake executable path: /usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/cmake-3.3.1-br2vzexkvky5litf32twjve4t4sescuq/bin/cmake

#######
# using gcc@4.7.1 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.7.1/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.7.1/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/apps/gnu/4.7.1/bin/gfortran" CACHE PATH "")

set(ENABLE_HDF5 ON CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/hdf5-1.8.16-tncq537skbm6gl3cfysngws6taadz2o7" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/conduit-github-3iz3pcoovcxdlv2sqcr4mveyqx4msenr" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/doxygen-1.8.10-n5kpsby5ccwqvkjxmuua72x35cddv6b5/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/python-2.7.8-g32izsb5gfir6ll5tsqa7hmcbvrpdblj/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/python-2.7.8-g32izsb5gfir6ll5tsqa7hmcbvrpdblj/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/uncrustify-0.61-v65xp2l6i423faxv4oa7amaqykzjaa7e/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/sparsehash-headers-2.0.2-ebexhuxsncs6fbkmof6djiv3d26rclik" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/boost-headers-1.58.0-zhkwg3db5a6xbucdytnjq3my52l5jlwu" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/lcov-1.11-tvagnsw74h4kwgos22odgtzjcbegi2xt/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.7.1/lcov-1.11-tvagnsw74h4kwgos22odgtzjcbegi2xt/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

#######
# MPI - added manually for now.
#######
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpif90" CACHE PATH "")
