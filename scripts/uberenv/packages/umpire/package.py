# Copyright 2013-2018 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import os


class Umpire(CMakePackage):
    """An application-focused API for memory management on NUMA & GPU
    architectures"""

    homepage = 'https://github.com/LLNL/Umpire'
    url='https://github.com/LLNL/Umpire/releases/download/v2.1.0/umpire-2.1.0.tar.gz'

    version('2.1.0', 'ec779f5a241c39a50d75eef1f248f61f')
    version('2.0.0', '4ba1c960b16552e0a85f28c244a9c383')
    version('1.0.0', '3435e9c9d48f2b05f838eb64fbf66d8e')
    version('0.3.2', '2b4138fb4d4272de8c34c073b4c6e88c')

    variant('cuda', default=False, description='Build with CUDA support')
    variant('fortran', default=False, description='Build C/Fortran API')
    variant('openmp', default=True, description='Build with OpenMP support')

    depends_on('cuda', when='+cuda')
    depends_on('cmake@3.3:', type='build')

    def cmake_args(self):
        spec = self.spec

        options = []

        if 'bgq' in os.getenv('SYS_TYPE', ""):
            options.extend(['-DENABLE_BENCHMARKS:BOOL=OFF',
                            '-DENABLE_TESTS:BOOL=OFF',
                            '-DENABLE_EXAMPLES:BOOL=OFF',
                            '-DENABLE_OPENMP:BOOL=OFF',
                            '-DBLT_CXX_FLAGS=-stdlib=libc++'])

        if '+cuda' in spec:
            options.extend([
                '-DENABLE_CUDA=On',
                '-DENABLE_DEVICE_CONST=On',
                '-DCUDA_TOOLKIT_ROOT_DIR=%s' % (spec['cuda'].prefix)])
        else:
            options.append('-DENABLE_CUDA=Off')

        if '+fortran' in spec:
            options.append('-DENABLE_FORTRAN=On')

        if '+openmp' in spec:
            options.append('-DENABLE_OPENMP=On')
        else:
            options.append('-DENABLE_OPENMP=Off')

        return options
