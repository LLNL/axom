# Copyright 2013-2018 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *


class Umpire(CMakePackage):
    """An application-focused API for memory management on NUMA & GPU
    architectures"""

    homepage = 'https://github.com/LLNL/Umpire'
    url='https://github.com/LLNL/Umpire/releases/download/v0.3.2/umpire-0.3.2.tar.gz'

    version('0.3.2', '2b4138fb4d4272de8c34c073b4c6e88c')

    variant('cuda', default=False, description='Build with CUDA support')
    variant('fortran', default=False, description='Build C/Fortran API')

    depends_on('cuda', when='+cuda')
    depends_on('cmake@3.3:', type='build')

    def cmake_args(self):
        spec = self.spec

        options = []

        if 'bgq' in env["SYS_TYPE"]:
            options.extend(['-DENABLE_BENCHMARKS:BOOL=OFF', 
                            '-DENABLE_TESTS:BOOL=OFF',
                            '-DENABLE_EXAMPLES:BOOL=OFF',
                            '-DENABLE_OPENMP:BOOL=OFF',
                            '-DBLT_CXX_FLAGS=-stdlib=libc++'])

        if '+cuda' in spec:
            options.extend([
                '-DENABLE_CUDA=On',
                '-DCUDA_TOOLKIT_ROOT_DIR=%s' % (spec['cuda'].prefix)])
        else:
            options.append('-DENABLE_CUDA=Off')

        if '+fortran' in spec:
            options.append('-DENABLE_FORTRAN=On')

        return options