##################################
# uberenv host-config
##################################
# x86_64-gcc@4.9.3
##################################

# cmake from uberenv
# cmake exectuable path: /home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/cmake-3.2.2-jp3rjtj6rx7upaun3ouujnhmkumzacfz/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/home/taylor16/gapps/gcc-4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/home/taylor16/gapps/gcc-4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/home/taylor16/gapps/gcc-4.9.3/bin/gfortran" CACHE PATH "")
set(GCOV_PATH "/home/taylor16/gapps/gcc-4.9.3/bin/gcov" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/conduit-github-t7wftaakgdvaxjbugnx7wv7z32a5e7ou" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/python-2.7.8-m7fqouw3e6rto2ygv7jjw6i5cofgc66c/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/python-2.7.8-m7fqouw3e6rto2ygv7jjw6i5cofgc66c/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/uncrustify-0.61-dgav63d2xyg7t3u4fjfvp5bjjxj5q245/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/sparsehash-headers-2.0.2-rkdejqec5amcvuorjml3qdcjtwazyvtk" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/boost-headers-1.58.0-byjp3xy3xzpwwratltskqsuyxqlk2vem" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/lcov-1.11-w5fedonfxcrbg2dycv5g4st4qdqh6pvv/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/home/taylor16/tpl/v2/spack/opt/spack/x86_64/gcc-4.9.3/lcov-1.11-w5fedonfxcrbg2dycv5g4st4qdqh6pvv/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

set(ENABLE_PYTHON_MODULE TRUE CACHE PATH "")
