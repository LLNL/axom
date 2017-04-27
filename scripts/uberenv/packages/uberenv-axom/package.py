from spack import *

import socket
import os
import platform
from os.path import join as pjoin

def cmake_cache_entry(name,value):
    return 'set(%s "%s" CACHE PATH "")\n\n' % (name,value)

class UberenvAxom(Package):
    """Spack Based Uberenv Build for Axom TPLs """

    # hash for dummy tarfile
    version('0.1', '8d378ef62dedc2df5db447b029b71200')
 
    homepage = "http://lc.llnl.gov/axom"

    # variants that allow us to winnow what TPLS we build
    # use ~cmake to skip cmake build and use whatever cmake is in the
    # users path 
    # (given the pain of building cmake on BGQ, this is really on for BGQ)
    variant('cmake',   default=True, description="Build cmake.")
    variant('devtools',default=True, description="Build development tools (such as sphinx, uncrustify, etc)")

    variant("python",  default=True, description="Build python")
    variant("lua",     default=True, description="Build lua")
    variant("boost",   default=True, description="Build boost headers")

    depends_on("conduit~shared+cmake",when="+cmake")
    depends_on("conduit~shared~cmake",when="~cmake")
    
    depends_on("hdf5~cxx~shared~fortran")

    # optional tpl builds
    depends_on("cmake@3.3.1",when="+cmake")

    depends_on("python",    when="+devtools")
    depends_on("doxygen",   when="+devtools")
    depends_on("uncrustify",when="+devtools")

    depends_on("python",   when="+python")
    depends_on("lua@5.1.5",when="+lua")

    # boost, header only
    depends_on("boost-headers", when="+boost")

    depends_on("py-sphinx", when="+devtools")
    depends_on("py-breathe",when="+devtools")
    depends_on("py-pyyaml", when="+devtools")
    depends_on("py-parsley",when="+devtools")

    if "darwin" in platform.system().lower():
        depends_on("mpich")

    if not "darwin" in platform.system().lower():
        depends_on("lcov",when="+devtools")

    # use dummy tarfile to avoid downloads
    def url_for_version(self, version):
        dummy_tar_path =  os.path.abspath(pjoin(os.path.dirname(__file__)))
        dummy_tar_path = pjoin(dummy_tar_path,"uberenv-axom.tar.gz")
        url      = "file://" + dummy_tar_path
        return url

    def install(self, spec, prefix):
        dest_dir     = env["SPACK_DEBUG_LOG_DIR"]
        
        c_compiler   = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]
        f_compiler   = None
        
        # see if we should enable fortran support
        if "SPACK_FC" in env.keys():
            # even if this is set, it may not exist
            # do one more sanity check
            if os.path.isfile(env["SPACK_FC"]):
                f_compiler  = env["SPACK_FC"]

        sys_type = spec.architecture
        # if on llnl systems, we can use the SYS_TYPE
        if env.has_key("SYS_TYPE"):
            sys_type = env["SYS_TYPE"]

        # cmake
        if "+cmake" in spec:
            cmake_exe = pjoin(spec['cmake'].prefix.bin,"cmake")
        else:
            cmake_exe = which("cmake")
            if cmake_exe is None:
                #error could not find cmake!
                crash()
            cmake_exe = cmake_exe.command
        
        host_cfg_fname = "%s-%s-%s.cmake" % (socket.gethostname().rstrip('1234567890'),sys_type,spec.compiler)
        host_cfg_fname = pjoin(dest_dir,host_cfg_fname)
        cfg = open(host_cfg_fname,"w")
        cfg.write("##################################\n")
        cfg.write("# uberenv host-config\n")
        cfg.write("#\n")
        cfg.write("# This is a generated file, edit at own risk.\n")
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
        cfg.write(cmake_cache_entry("CMAKE_C_COMPILER",c_compiler))

        cfg.write("# cpp compiler used by spack\n")
        cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER",cpp_compiler))
        
        cfg.write("# fortran compiler used by spack\n")
        if not f_compiler is None:
            cfg.write(cmake_cache_entry("ENABLE_FORTRAN","ON"))
            cfg.write(cmake_cache_entry("CMAKE_Fortran_COMPILER",f_compiler))
        else:
            cfg.write("# no fortran compiler\n\n")
            cfg.write(cmake_cache_entry("ENABLE_FORTRAN","OFF"))

        # try to find the common prefix of the TPL directory, including the compiler 
        # Fallback is to set tpl_root to full install path for uberenv-axom, which will not be in other paths
        compiler_str = str(spec.compiler).replace('@','-')
        prefix_paths = prefix.split( compiler_str )

        if len(prefix_paths) == 2:
            tpl_root = pjoin( prefix_paths[0], compiler_str )
            cfg.write("# Root directory for generated TPLs\n")
            cfg.write(cmake_cache_entry("TPL_ROOT",tpl_root))
        else:
            tpl_root = prefix  

        hdf5_dir = spec["hdf5"].prefix.replace(tpl_root, "${TPL_ROOT}")
        cfg.write("# hdf5 from uberenv\n")
        cfg.write(cmake_cache_entry("HDF5_DIR",hdf5_dir))

        conduit_dir = spec["conduit"].prefix.replace(tpl_root, "${TPL_ROOT}")
        cfg.write("# conduit from uberenv\n")
        cfg.write(cmake_cache_entry("CONDUIT_DIR",conduit_dir))

        if "boost-headers" in spec:
            boost_headers_dir = spec["boost-headers"].prefix.replace(tpl_root, "${TPL_ROOT}")
            cfg.write("# boost headers from uberenv\n")
            cfg.write(cmake_cache_entry("BOOST_DIR",boost_headers_dir))
        else:
            cfg.write("# boost headers not installed by uberenv\n")

        # optional tpls

        if "python" in spec or "devtools" in spec:
            python_exe = pjoin(spec['python'].prefix.bin,"python").replace(tpl_root, "${TPL_ROOT}")
            cfg.write("# python from uberenv\n")
            cfg.write(cmake_cache_entry("PYTHON_EXECUTABLE",python_exe))
        else:
            cfg.write("# python not build by uberenv\n\n")

        if "lua" in spec:
            lua_dir = spec["lua"].prefix.replace(tpl_root, "${TPL_ROOT}")
            cfg.write("# lua from uberenv\n")
            cfg.write(cmake_cache_entry("LUA_DIR",lua_dir))
        else:
            cfg.write("# lua not build by uberenv\n\n")

        # optional tpls (dev tools)

        if "doxygen" in spec:
            doxygen_exe = pjoin(spec['doxygen'].prefix.bin,"doxygen").replace(tpl_root, "${TPL_ROOT}")
            cfg.write("# doxygen from uberenv\n")
            cfg.write(cmake_cache_entry("DOXYGEN_EXECUTABLE", doxygen_exe))
        else:
            cfg.write("# doxygen not built by uberenv\n\n")

        if "py-sphinx" in spec:
            sphinx_build_exe = pjoin(spec['python'].prefix.bin,"sphinx-build").replace(tpl_root, "${TPL_ROOT}")
            cfg.write("# sphinx from uberenv\n")
            cfg.write(cmake_cache_entry("SPHINX_EXECUTABLE", sphinx_build_exe))
        else:
            cfg.write("# sphinx not built by uberenv\n\n")

        if "uncrustify" in spec:
            uncrustify_exe   = pjoin(spec['uncrustify'].prefix.bin,"uncrustify").replace(tpl_root, "${TPL_ROOT}")
            cfg.write("# uncrustify from uberenv\n")
            cfg.write(cmake_cache_entry("UNCRUSTIFY_EXECUTABLE", uncrustify_exe))
        else:
            cfg.write("# uncrustify not built by uberenv\n\n")

        cfg.write("# lcov and genhtml from uberenv\n")
        if "lcov" in spec:
            lcov_path = spec['lcov'].prefix.replace(tpl_root, "${TPL_ROOT}")
            cfg.write(cmake_cache_entry("LCOV_PATH", pjoin(lcov_path,"usr","bin","lcov")))
            cfg.write(cmake_cache_entry("GENHTML_PATH",pjoin(lcov_path,"usr","bin","genhtml")))
        else:
            cfg.write("# lcov and genhtml not built by uberenv\n\n")

        #MPI SUPPORT (when mpich is enabled)
        if "mpich" in spec:
            mpiexec  = pjoin(spec['mpich'].prefix.bin,"mpiexec")
            mpicc    = pjoin(spec['mpich'].prefix.bin,"mpicc")
            mpif90   = pjoin(spec['mpich'].prefix.bin,"mpif90")
            cfg.write("# MPI Support\n")
            cfg.write(cmake_cache_entry("ENABLE_MPI","ON"))
            cfg.write(cmake_cache_entry("MPI_C_COMPILER",mpicc))
            # we use `mpicc` as `MPI_CXX_COMPILER` b/c we don't want to introduce 
            # linking deps to the MPI C++ libs (we aren't using C++ features of MPI)
            cfg.write(cmake_cache_entry("MPI_CXX_COMPILER",mpicc))
            cfg.write(cmake_cache_entry("MPI_Fortran_COMPILER",mpif90))
            cfg.write(cmake_cache_entry("MPIEXEC",mpiexec))


        # Note (KW 3/2016) -- per ATK-659, we are temporarily disabling CXX11 for default configurations on intel builds 
        if "intel" in spec.compiler.name:
            cfg.write("# Temporarily disable CXX11 on intel builds until we resolve issue ATK-619\n")
            cfg.write(cmake_cache_entry("BLT_CXX_STD","c++98"))

        cfg.write("##################################\n")
        cfg.write("# end uberenv host-config\n")
        cfg.write("##################################\n")
        cfg.write("\n")
        cfg.close()
        mkdirp(prefix)
        install(host_cfg_fname,prefix)
        print "[result host-config file: %s]" % host_cfg_fname



