##################################
# uberenv host-config
##################################
# x86_64-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/cmake-3.3.1-jkk5puodwd24sh2egavyrp2ony2z3wtn/bin/cmake

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

# conduit from uberenv
set(CONDUIT_DIR "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/conduit-github-4nhhkakacu5nismusjvmqxxyhipmgdmy" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/doxygen-1.8.10-5bivwj24gorejk2bvq2ltvfkw5aty2mh/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/python-2.7.8-dqhz7ufkpnajumalaa4damywigthzlvi/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/python-2.7.8-dqhz7ufkpnajumalaa4damywigthzlvi/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/uncrustify-0.61-sdasmocwkfmf6xupjpoo73rzdddf2nl6/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/sparsehash-headers-2.0.2-rmpk4wmjkozch4i54wa7l5zjfyru7yr3" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/boost-headers-1.58.0-6jjq5njm5itr67imya3p4fcnmiaj5nab" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/lcov-1.11-h3ijl43wxo53yrr7j7gsnmpuc4wc5nrn/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/home/taylor16/tpl/v4/spack/opt/spack/x86_64/gcc-4.9.3/lcov-1.11-h3ijl43wxo53yrr7j7gsnmpuc4wc5nrn/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

