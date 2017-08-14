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
    homepage = "http://software.llnl.gov/conduit/"
    url      = "https://github.com/LLNL/conduit/archive/v0.2.1.tar.gz"

    version('0.2.1', 'cd2b42c76f70ac3546582b6da77c6028')
    version('0.2.0', 'd595573dedf55514c11d7391092fd760')

    variant('cmake',  default=True, description="Build cmake.")
    variant('shared', default=True, description="Build shared libraries.")
    
    depends_on("cmake@3.8.2",when="+cmake")
    depends_on("hdf5~shared~fortran")
    
    if "darwin" in platform.system().lower():
        depends_on("mpich")


    def install(self, spec, prefix):
        with working_dir('spack-build', create=True):
            hdf5_dir = spec['hdf5'].prefix
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

            cmake_args.append("-DHDF5_DIR=%s" % hdf5_dir);
            # turn off docs, so we don't have problems 
            # w/ system sphinx installs
            cmake_args.append("-DENABLE_DOCS=OFF") 
            
            # see if we should enable fortran support
            f_compiler   = None
            if "SPACK_FC" in env.keys():
                if os.path.isfile(env["SPACK_FC"]):
                    f_compiler = env["SPACK_FC"]
                    cmake_args.append("-DENABLE_FORTRAN=ON")


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



  
