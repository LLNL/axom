##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-intel@16.0.109
##################################

# cmake from uberenv
# cmake executable path: /usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/cmake-3.3.1-b7plebpo7hhjv5owhziowxshxvu6ebe3/bin/cmake

#######
# using intel@16.0.109 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/tools/ic-16.0.109/bin/icc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/tools/ic-16.0.109/bin/icpc" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/local/tools/ic-16.0.109/bin/ifort" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/conduit-github-hodnxg7yb62xzouc265rf4pjhhhppnrk" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/doxygen-1.8.10-i6mt67xh7b7qjfg76gnbjurhgymfq7va/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/python-2.7.8-ssya74k3fc7yjdasoqcpkqgt7xukukk3/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/python-2.7.8-ssya74k3fc7yjdasoqcpkqgt7xukukk3/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/uncrustify-0.61-7v4io4dwzqwtplzuykvomoxpibd3cnbl/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/sparsehash-headers-2.0.2-cv3sza7qef4mtf4lthelmarncfaz7ygl" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/boost-headers-1.58.0-vaxfagykvpueronnzybsyr63rhcmldzj" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/lcov-1.11-nrsx4h6lyz4ajjrqi45dphsnbppj3da3/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/stable/spack/opt/spack/chaos_5_x86_64_ib/intel-16.0.109/lcov-1.11-nrsx4h6lyz4ajjrqi45dphsnbppj3da3/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

#######
# MPI 
#######
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/mvapich2-intel-2.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/mvapich2-intel-2.0/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER  "/usr/local/tools/mvapich2-intel-2.0/bin/mpif90" CACHE PATH "")
