# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Raja(CMakePackage):
    """RAJA Parallel Framework."""

    homepage = "http://software.llnl.gov/RAJA/"
    git      = "https://github.com/LLNL/RAJA.git"

    version('develop', branch='develop', submodules='True')
    version('master',  branch='master',  submodules='True')
    version('0.11.0', tag='v0.11.0', submodules="True")
    version('0.10.1', tag='v0.10.1', submodules="True")
    version('0.9.0', tag='v0.9.0', submodules="True")
    version('0.8.0', tag='v0.8.0', submodules="True")
    version('0.7.0', tag='v0.7.0', submodules="True")
    version('0.6.0', tag='v0.6.0', submodules="True")
    version('0.5.3', tag='v0.5.3', submodules="True")
    version('0.5.2', tag='v0.5.2', submodules="True")
    version('0.5.1', tag='v0.5.1', submodules="True")
    version('0.5.0', tag='v0.5.0', submodules="True")
    version('0.4.1', tag='v0.4.1', submodules="True")
    version('0.4.0', tag='v0.4.0', submodules="True")

    variant('cuda', default=False, description='Build with CUDA backend')
    variant('openmp', default=True, description='Build OpenMP backend')

    depends_on('cuda', when='+cuda')

    depends_on('cmake@3.8:', type='build')
    depends_on('cmake@3.9:', when='+cuda', type='build')

    def _get_sys_type(self, spec):
        sys_type = spec.architecture
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        return sys_type

    def cmake_args(self):
        spec = self.spec

        sys_type = self._get_sys_type(spec)
        on_blueos = 'blueos' in sys_type
        on_blueos_p9 = on_blueos and 'p9' in sys_type

        options = []
        options.append('-DENABLE_OPENMP={0}'.format(
            'On' if '+openmp' in spec else 'Off'))

        if '+cuda' in spec:
            if on_blueos_p9:
                options.extend(['-DCMAKE_CUDA_FLAGS:STRING=-arch sm_70'])
            elif on_blueos:
                options.extend(['-DCMAKE_CUDA_FLAGS:STRING=-arch sm_60'])

            options.extend([
                '-DENABLE_CUDA=On',
                '-DCUDA_TOOLKIT_ROOT_DIR=%s' % (spec['cuda'].prefix)])

        # Work around spack adding -march=ppc64le to SPACK_TARGET_ARGS which is used by the spack compiler wrapper
        if on_blueos:
            options.extend(['-DENABLE_TESTS=OFF'])

        return options
