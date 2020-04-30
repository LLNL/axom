# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *


class Umpire(CMakePackage):
    """An application-focused API for memory management on NUMA & GPU
    architectures"""

    homepage = 'https://github.com/LLNL/Umpire'
    git      = 'https://github.com/LLNL/Umpire.git'

    version('develop', branch='develop', submodules='True')
    version('master', branch='master', submodules='True')
    version('2.1.0', tag='v2.1.0', submodules='True')
    version('2.0.0', tag='v2.0.0', submodules='True')
    version('1.1.0', tag='v1.1.0', submodules='True')
    version('1.0.1', tag='v1.0.1', submodules='True')
    version('1.0.0', tag='v1.0.0', submodules='True')
    version('0.3.5', tag='v0.3.5', submodules='True')
    version('0.3.3', tag='v0.3.3', submodules='True')
    version('0.3.2', tag='v0.3.2', submodules='True')
    version('0.3.1', tag='v0.3.1', submodules='True')
    version('0.3.0', tag='v0.3.0', submodules='True')
    version('0.2.4', tag='v0.2.4', submodules='True')
    version('0.2.3', tag='v0.2.3', submodules='True')
    version('0.2.2', tag='v0.2.2', submodules='True')
    version('0.2.1', tag='v0.2.1', submodules='True')
    version('0.2.0', tag='v0.2.0', submodules='True')
    version('0.1.4', tag='v0.1.4', submodules='True')
    version('0.1.3', tag='v0.1.3', submodules='True')

    variant('cuda', default=False, description='Build with CUDA support')
    variant('fortran', default=False, description='Build C/Fortran API')
    variant('numa', default=False, description='Enable NUMA support')
    variant('openmp', default=False, description='Build with OpenMP support')

    depends_on('cuda', when='+cuda')
    depends_on('cmake@3.8:', type='build')
    depends_on('cmake@3.9:', when='+cuda', type='build')

    conflicts('+numa', when='@:0.3.2')

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

        if '+cuda' in spec:
            if on_blueos_p9:
                options.extend(['-DCMAKE_CUDA_FLAGS:STRING=-arch sm_70'])
            elif on_blueos:
                options.extend(['-DCMAKE_CUDA_FLAGS:STRING=-arch sm_60'])

            options.extend([
                '-DENABLE_CUDA=On',
                '-DENABLE_DEVICE_CONST=On',
                '-DCUDA_TOOLKIT_ROOT_DIR=%s' % (spec['cuda'].prefix)])
        else:
            options.append('-DENABLE_CUDA=Off')

        if '+fortran' in spec:
            options.append('-DENABLE_FORTRAN=On')

        if '+numa' in spec:
            options.append('-DENABLE_NUMA=On')

        if '+openmp' in spec:
            options.append('-DENABLE_OPENMP=On')

        return options
