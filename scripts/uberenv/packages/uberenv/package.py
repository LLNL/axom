from spack import *

import socket
from os.path import join as pjoin

class Uberenv(Package):
    """Spack Based Uber Build for Toolkit TPLs """

    # TODO: what do we need for DIY mode?
    
    homepage = "http://sphinx-doc.org/"
    url      = "https://pypi.python.org/packages/source/S/Sphinx/Sphinx-1.3.1.tar.gz#md5=8786a194acf9673464c5455b11fd4332"

    version('1.3.1', '8786a194acf9673464c5455b11fd4332')
    
    # both of theses are custom
    depends_on("python")
    depends_on("py-sphinx")
    depends_on("py-breathe")
    
    # this was pushed to develop, but not yet in the diy branch
    depends_on("uncrustify")

    # boost, header only
    depends_on("boost-headers")
        
    def install(self, spec, prefix):
        #mkdirp(prefix)
        dest_dir = env["SPACK_DEBUG_LOG_DIR"]
        sphinx_build_exe = pjoin(spec['python'].prefix.bin,"sphinx-build")
        uncrustify_exe   = pjoin(spec['uncrustify'].prefix.bin,"uncrustify")
        # TODO: better name (use sys-type and compiler name ?)
        cfg = open(pjoin(dest_dir,"%s.cmake" % socket.gethostname()),"w")
        cfg.write("#######\n")
        cfg.write("# uberenv host-config for asctoolkit\n")
        cfg.write("#######\n")
        cfg.write("# sphinx from uberenv\n")
        cfg.write('set(SPHINX_EXECUTABLE "%s" CACHE PATH "")\n\n' % sphinx_build_exe)
        cfg.write("# uncrustify from uberenv\n")
        cfg.write('set(UNCRUSTIFY_EXECUTABLE "%s" CACHE PATH "")\n\n' % uncrustify_exe)
        cfg.write("# boost headers from uberenv\n")
        cfg.write('set(BOOST_ROOT "%s" CACHE PATH "")\n\n' % spec['boost-headers'].prefix)
        cfg.write("\n")
        cfg.close()        
        
        
        