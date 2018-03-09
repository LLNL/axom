###############################################################################
# Copyright (c) 2014-2015, Lawrence Livermore National Security, LLC.
# 
# Produced at the Lawrence Livermore National Laboratory
# 
# LLNL-CODE-666778
# 
# All rights reserved.
# 
# This file is part of Conduit. 
# 
# For details, see https://lc.llnl.gov/conduit/.
# 
# Please also read conduit/LICENSE
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the disclaimer below.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
# 
###############################################################################
from spack import *

import os
import sys
import platform
from os.path import join as pjoin

class Conduit(Package):
    """Conduit is an open source project from Lawrence Livermore National
    Laboratory that provides an intuitive model for describing hierarchical
    scientific data in C++, C, Fortran, and Python. It is used for data
    coupling between packages in-core, serialization, and I/O tasks."""

    homepage = "http://software.llnl.gov/conduit"
    url = "https://github.com/LLNL/conduit/releases/download/v0.3.0/conduit-v0.3.0-src-with-blt.tar.gz"

    version('0.3.1', 'b98d1476199a46bde197220cd9cde042')
    version('0.3.0', '6396f1d1ca16594d7c66d4535d4f898e')
    # note: checksums on github automatic release source tars changed ~9/17
    version('0.2.1', 'ed7358af3463ba03f07eddd6a6e626ff')
    version('0.2.0', 'a7b398d493fd71b881a217993a9a29d4')

    variant('cmake',  default=True, description="Build cmake.")
    variant('shared', default=True, description="Build shared libraries.")
    variant("hdf5",   default=True, description="Build hdf5")
    
    # cmake 3.8.2 or newer
    depends_on("cmake@3.8.2:",when="+cmake")
    # hdf5 1.8.x
    depends_on("hdf5@1.8.16:1.8.999~shared~fortran", when="+hdf5")
    
    if "darwin" in platform.system().lower():
        depends_on("mpich")

    def url_for_version(self, version):
        """
        Provide proper url
        """
        v = str(version)
        if v == "0.2.0":
            return "https://github.com/LLNL/conduit/archive/v0.2.0.tar.gz"
        elif v == "0.2.1":
            return "https://github.com/LLNL/conduit/archive/v0.2.1.tar.gz"
        else:
            # starting with v 0.3.0, conduit uses BLT
            # (https://github.com/llnl/blt) as a submodule, since github does
            # not automatically package source from submodules, conduit
            # provides a custom src tarball
            return "https://github.com/LLNL/conduit/releases/download/v{0}/conduit-v{1}-src-with-blt.tar.gz".format(v, v)
        return url

    def install(self, spec, prefix):
        with working_dir('spack-build', create=True):
            cmake_args = ["../src"]
            # if we have a static build, we need to avoid
            # any of cmake's default settings related to 
            # rpath
            if "+shared" in spec:
                cmake_args.extend(std_cmake_args)
            else:
                for arg in std_cmake_args:
                    if arg.count("RPATH") == 0:
                        cmake_args.append(arg)

            # Add HDF5 directory and fix HDF5 library paths on BGQ
            # to deal with issues with cmake's FindHDF5 on that platform
            if "+hdf5" in spec:
                hdf5_dir = spec['hdf5'].prefix
                cmake_args.append("-DHDF5_DIR=%s" % hdf5_dir)
                if 'bgqos_0' in os.getenv('SYS_TYPE', ""):
                    cmake_args.append("-DHDF5_C_LIBRARY_m:STRING=-lm")
                    cmake_args.append("-DHDF5_C_LIBRARY_dl:STRING=-ldl")
                
            # turn off docs, so we don't have problems 
            # w/ system sphinx installs
            cmake_args.append("-DENABLE_DOCS=OFF") 
            
            # see if we should enable fortran support
            f_compiler   = None
            if "SPACK_FC" in env.keys():
                if os.path.isfile(env["SPACK_FC"]):
                    f_compiler = env["SPACK_FC"]
                    cmake_args.append("-DENABLE_FORTRAN=ON")
                    
                    # Fix missing std linker flag in xlc compiler
                    # TODO: Fix to remove hard-coded compiler path
                    if 'blueos_3' in os.getenv('SYS_TYPE', "") and 'xl@coral' in os.getenv('SPACK_COMPILER_SPEC', ""):
                        cmake_args.append("-DBLT_FORTRAN_FLAGS:STRING=-WF,-C! -Wl,-lstdc++ -Wl,/usr/tce/packages/xl/xl-beta-2017.10.13/xlC/13.1.6/lib/libibmc++.so")

            if "+shared" in spec:
                cmake_args.append("-DBUILD_SHARED_LIBS=ON")
            else:
                cmake_args.append("-DBUILD_SHARED_LIBS=OFF")

            if "mpich" in spec:
                mpiexec  = pjoin(spec['mpich'].prefix.bin,"mpiexec")
                mpicc    = pjoin(spec['mpich'].prefix.bin,"mpicc")
                cmake_args.append("-DENABLE_MPI=ON")
                cmake_args.append("-DMPI_C_COMPILER:PATH=%s"   % mpicc)
                cmake_args.append("-DMPI_CXX_COMPILER:PATH=%s" % mpicc)

                cmake_args.append("-DMPIEXEC:PATH=%s" % mpiexec)
                
                if not f_compiler is None:
                    mpif90   = pjoin(spec['mpich'].prefix.bin,"mpif90")
                    cmake_args.append("-DMPI_Fortran_COMPILER:PATH=%s" % mpif90)
                
            cmake(*cmake_args)
            make(parallel=False)
            make("install",parallel=False)



  
