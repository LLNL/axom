from spack import *

import socket
import os
import platform
from os.path import join as pjoin

class ToolkitCpp11MinimalTpls(Package):
    """Spack package for minimal set of TPLs for toolkit C++11 build, no docs or python """

    # hash for dummy tarfile
    version('0.1', '37b303be8b9da35e9369bb8599237256')
 
    homepage = "http://lc.llnl.gov/toolkit"
    # standard spack packages
    depends_on("hdf5~shared~zlib")

    # custom spack packages
    depends_on("cmake~ncurses~openssl@3.3.1")
    depends_on("conduit~shared")

    # use dummy tarfile to avoid downloads
    def url_for_version(self, version):
        dummy_tar_path =  os.path.abspath(pjoin(os.path.dirname(__file__)))
        dummy_tar_path = pjoin(dummy_tar_path,"toolkit-cpp11-minimal-tpls-0.1.tar.gz")
        url      = "file://" + dummy_tar_path
        return url

    def install(self, spec, prefix):
        #mkdirp(prefix)
        dest_dir     = env["SPACK_DEBUG_LOG_DIR"]
        c_compiler   = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]
        f_compiler = None
        if "SPACK_FC" in env.keys():
            f_compiler   = env["SPACK_FC"]
        sys_type     = spec.architecture
        if "SYS_TYPE" in env.keys():
            sys_type = env["SYS_TYPE"]
        # conduit
        conduit_dir      = spec['conduit'].prefix
        cmake_exe        = pjoin(spec['cmake'].prefix.bin,"cmake")

        host_cfg_fname = "%s-%s-%s.cmake" % (socket.gethostname().rstrip('1234567890'),sys_type,spec.compiler)
        host_cfg_fname = pjoin(dest_dir,host_cfg_fname)
        cfg = open(host_cfg_fname,"w")
        cfg.write("##################################\n")
        cfg.write("# uberenv host-config\n")
        cfg.write("##################################\n")
        cfg.write("# %s-%s\n" % (sys_type,spec.compiler))
        cfg.write("##################################\n\n")
        # show path to cmake for reference
        cfg.write("# cmake from uberenv\n")
        cfg.write("# cmake executable path: %s\n\n" % cmake_exe)

        # compiler settings
        cfg.write("#######\n")
        cfg.write("# using %s compiler spec\n" % spec.compiler)
        cfg.write("#######\n\n")
        cfg.write("# c compiler used by spack\n")
        cfg.write('set(CMAKE_C_COMPILER "%s" CACHE PATH "")\n\n' % c_compiler)
        cfg.write("# cpp compiler used by spack\n")
        cfg.write('set(CMAKE_CXX_COMPILER "%s" CACHE PATH "")\n\n' % cpp_compiler)
        cfg.write("# fortran compiler used by spack\n")
        if not f_compiler is None:
            cfg.write('set(ENABLE_FORTRAN ON CACHE PATH "")\n\n')
            cfg.write('set(CMAKE_Fortran_COMPILER  "%s" CACHE PATH "")\n\n' % f_compiler)
        else:
            cfg.write("# no fortran compiler\n\n")
            cfg.write('set(ENABLE_FORTRAN OFF CACHE PATH "")\n\n')

        cfg.write('set(ENABLE_HDF5 ON CACHE PATH "")\n\n')
        cfg.write("# hdf5 from uberenv\n")
        cfg.write('set(HDF5_DIR "%s" CACHE PATH "")\n\n' % spec['hdf5'].prefix)

        cfg.write("# conduit from uberenv\n")
        cfg.write('set(CONDUIT_DIR "%s" CACHE PATH "")\n\n' % conduit_dir)

        cfg.write("##################################\n")
        cfg.write("# end uberenv host-config\n")
        cfg.write("##################################\n")
        cfg.write("\n")
        cfg.close()
        mkdirp(prefix)
        install(host_cfg_fname,prefix)
        print "[result host-config file: %s]" % host_cfg_fname
