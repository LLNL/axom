# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Raja(CMakePackage):
    """RAJA Parallel Framework."""

    homepage = "http://software.llnl.gov/RAJA/"
    url='https://github.com/LLNL/RAJA/releases/download/v0.7.0/RAJA-0.7.0.tar.gz'

    version('0.7.0', '18b9654d75ee2f9278b3a8f19f9989fa')

    variant('cuda', default=False, description='Build with CUDA backend')
    variant('openmp', default=True, description='Build OpenMP backend')

    depends_on('cuda', when='+cuda')

    depends_on('cmake@3.3:', type='build')

    def cmake_args(self):
        spec = self.spec

        options = []

        if 'bgq' in env["SYS_TYPE"]:
            options.extend(['-DBLT_CXX_FLAGS=-stdlib=libc++'])

        if '+openmp' in spec:
            options.extend([
                '-DENABLE_OPENMP=On'])
        else:
            options.extend([
                '-DENABLE_OPENMP=Off'])

        if '+cuda' in spec:
            options.extend([
                '-DENABLE_CUDA=On',
                '-DCUDA_TOOLKIT_ROOT_DIR=%s' % (spec['cuda'].prefix)])

        return options