from spack import *

import socket
import os
import platform
from os.path import join as pjoin

class UberenvAsctoolkit(Package):
    """Spack Based Uberenv Build for Toolkit TPLs """

    # hash for dummy tarfile
    version('0.1', '8d378ef62dedc2df5db447b029b71200')
    
    # all of theses are custom
    depends_on("cmake")
    depends_on("doxygen")
    depends_on("python")
    depends_on("conduit")
    depends_on("py-sphinx")
    depends_on("py-breathe")
    depends_on("py-pyyaml")
    depends_on("py-parsley")
    depends_on("py-cogapp")

    if not "darwin" in platform.system().lower():
        depends_on("lcov")

    # boost, header only
    depends_on("boost-headers")
    depends_on("sparsehash-headers")

    # this was pushed to develop, but not yet in the diy branch
    depends_on("uncrustify")

    # use dummy tarfile to avoid downloads
    def url_for_version(self, version):
        dummy_tar_path =  os.path.abspath(pjoin(os.path.split(__file__)[0]))
        dummy_tar_path = pjoin(dummy_tar_path,"uberenv-asctoolkit.tar.gz")
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
        doxygen_exe      = pjoin(spec['doxygen'].prefix.bin,"doxygen")
        python_exe       = pjoin(spec['python'].prefix.bin,"python")
        sphinx_build_exe = pjoin(spec['python'].prefix.bin,"sphinx-build")
        uncrustify_exe   = pjoin(spec['uncrustify'].prefix.bin,"uncrustify")

        host_cfg_fname = "%s-%s-%s.cmake" % (socket.gethostname(),sys_type,spec.compiler)
        host_cfg_fname = pjoin(dest_dir,host_cfg_fname)
        cfg = open(host_cfg_fname,"w")
        cfg.write("##################################\n")
        cfg.write("# uberenv host-config\n")
        cfg.write("##################################\n")
        cfg.write("# %s-%s\n" % (sys_type,spec.compiler))
        cfg.write("##################################\n\n")
        # show path to cmake for reference
        cfg.write("# cmake from uberenv\n")
        cfg.write("# cmake exectuable path: %s\n\n" % cmake_exe)

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

        cfg.write("# conduit from uberenv\n")
        cfg.write('set(CONDUIT_DIR "%s" CACHE PATH "")\n\n' % conduit_dir)

        cfg.write("# doxygen from uberenv\n")
        cfg.write('set(DOXYGEN_EXECUTABLE "%s" CACHE PATH "")\n\n' % doxygen_exe)

        cfg.write("# python from uberenv\n")
        cfg.write('set(PYTHON_EXECUTABLE "%s" CACHE PATH "")\n\n' % python_exe)
        
        cfg.write("# sphinx from uberenv\n")
        cfg.write('set(SPHINX_EXECUTABLE "%s" CACHE PATH "")\n\n' % sphinx_build_exe)

        cfg.write("# uncrustify from uberenv\n")
        cfg.write('set(UNCRUSTIFY_EXECUTABLE "%s" CACHE PATH "")\n\n' % uncrustify_exe)

        cfg.write("# sparsehash headers from uberenv\n")
        cfg.write('set(SPARSEHASH_DIR "%s" CACHE PATH "")\n\n' % spec['sparsehash-headers'].prefix)

        cfg.write("# boost headers from uberenv\n")
        cfg.write('set(ENABLE_BOOST ON CACHE PATH "")\n')
        cfg.write('set(BOOST_ROOT "%s" CACHE PATH "")\n\n' % spec['boost-headers'].prefix)

        cfg.write("# lcov and genhtml from uberenv\n")
        if "lcov" in spec:
            cfg.write('set(LCOV_PATH "%s/usr/bin/lcov" CACHE PATH "")\n\n' % spec['lcov'].prefix)
            cfg.write('set(GENHTML_PATH "%s/usr/bin/genhtml" CACHE PATH "")\n\n' % spec['lcov'].prefix)
        else:
            cfg.write("# lcov and genhtml not built by uberenv\n\n")

        cfg.write("##################################\n")
        cfg.write("# end uberenv host-config\n")
        cfg.write("##################################\n")
        cfg.write("\n")
        cfg.close()
        mkdirp(prefix)
        install(host_cfg_fname,prefix)
        print "[result host-config file: %s]" % host_cfg_fname




