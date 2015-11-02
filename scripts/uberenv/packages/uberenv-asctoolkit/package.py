from spack import *

import socket
import os
import platform
from os.path import join as pjoin

class UberenvAsctoolkit(Package):
    """Spack Based Uberenv Build for Toolkit TPLs """

    # hash for dummy tarfile
    version('0.1', '8d378ef62dedc2df5db447b029b71200')
    
    # all of these theses are custom
    depends_on("cmake")
    depends_on("python")
    depends_on("conduit")
    depends_on("py-sphinx")
    depends_on("py-breathe")
    depends_on("py-pyyaml")
    depends_on("py-parsley")

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
        dest_dir = env["SPACK_DEBUG_LOG_DIR"]
        cmake_exe        = pjoin(spec['cmake'].prefix.bin,"cmake")
        conduit_dir      = spec['conduit'].prefix
        python_exe       = pjoin(spec['python'].prefix.bin,"python")
        sphinx_build_exe = pjoin(spec['python'].prefix.bin,"sphinx-build")
        uncrustify_exe   = pjoin(spec['uncrustify'].prefix.bin,"uncrustify")
        # TODO: better name (use sys-type and compiler name ?)
        cfg = open(pjoin(dest_dir,"%s.cmake" % socket.gethostname()),"w")
        cfg.write("#######\n")
        cfg.write("# uberenv host-config for asctoolkit\n")
        cfg.write("#######\n")
        cfg.write("# cmake from uberenv\n")
        cfg.write("# cmake exectuable path: %s\n\n" % cmake_exe)

        cfg.write("# conduit from uberenv\n")
        cfg.write('set(CONDUIT_DIR "%s" CACHE PATH "")\n\n' % conduit_dir)

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

        cfg.write("\n")
        cfg.close()
        
        
        
