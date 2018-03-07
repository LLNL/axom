from spack import *

import socket
import os
import platform
from os.path import join as pjoin

def cmake_cache_entry(name,value):
    """Generate a string for a cmake cache variable""" 
    return 'set(%s "%s" CACHE PATH "")\n\n' % (name,value)

def cmake_cache_option(name,boolean_value):
    """Generate a string for a cmake configuration option""" 
    
    value = "ON" if boolean_value else "OFF"
    return 'set(%s %s CACHE BOOL "")\n\n' % (name,value)

def get_spec_path(spec, package_name, path_replacements = {}, use_bin = False) :
    """Extracts the prefix path for the given spack package
       path_replacements is a dictionary with string replacements for the path.
    """

    if not use_bin:
        path = spec[package_name].prefix
    else:
        path = spec[package_name].prefix.bin
        
    for key in path_replacements:
        path = path.replace(key,path_replacements[key])
        
    return path


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

    variant("mfem",   default=True, description="Build mfem")

    variant("hdf5",   default=True, description="Build hdf5")

    variant("scr", default=False, description="Build SCR")

    depends_on("hdf5~cxx~shared~fortran", when="+hdf5")

    depends_on("conduit~shared") # needed so we can optionally use ^conduit@master

    depends_on("conduit~shared+cmake+hdf5",when="+cmake+hdf5")
    depends_on("conduit~shared~cmake+hdf5",when="~cmake+hdf5")
    depends_on("conduit~shared+cmake~hdf5",when="+cmake~hdf5")
    depends_on("conduit~shared~cmake~hdf5",when="~cmake~hdf5")
    
    depends_on("scr", when="+scr")

    # builds serial version of mfem that does not depend on Sidre
    depends_on("mfem~mpi~gzstream",   when="+mfem")

    # optional tpl builds
    depends_on("cmake@3.9.6",when="+cmake")

    depends_on("python",    when="+devtools")
    depends_on("doxygen",   when="+devtools")
    depends_on("uncrustify",when="+devtools")

    depends_on("python",   when="+python")

    depends_on("py-sphinx", when="+devtools")
    #depends_on("py-breathe",when="+devtools")
    depends_on("py-shroud", when="+devtools")

    if "darwin" in platform.system().lower():
        depends_on("mpich")

    #if not "darwin" in platform.system().lower():
    #    depends_on("lcov",when="+devtools")

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
            cfg.write(cmake_cache_option("ENABLE_FORTRAN",True))
            cfg.write(cmake_cache_entry("CMAKE_Fortran_COMPILER",f_compiler))
        else:
            cfg.write("# no fortran compiler\n\n")
            cfg.write(cmake_cache_option("ENABLE_FORTRAN",False))

        # Try to find the common prefix of the TPL directory, including the compiler 
        # If found, we will use this in the TPL paths
        compiler_str = str(spec.compiler).replace('@','-')
        prefix_paths = prefix.split( compiler_str )
        path_replacements = {}

        if len(prefix_paths) == 2:
            tpl_root = pjoin( prefix_paths[0], compiler_str )
            path_replacements[tpl_root] = "${TPL_ROOT}"
            cfg.write("# Root directory for generated TPLs\n")
            cfg.write(cmake_cache_entry("TPL_ROOT",tpl_root))

        if "hdf5" in spec:
            hdf5_dir = get_spec_path(spec, "hdf5", path_replacements)
            cfg.write("# hdf5 from uberenv\n")
            cfg.write(cmake_cache_entry("HDF5_DIR",hdf5_dir))
        else:
            cfg.write("# hdf5 not built by uberenv\n\n")

        if "scr" in spec:
            scr_dir = get_spec_path(spec, "scr", path_replacements)
            cfg.write("# scr from uberenv\n")
            cfg.write(cmake_cache_entry("SCR_DIR",scr_dir))
        else:
            cfg.write("# scr not built by uberenv\n\n")

        conduit_dir = get_spec_path(spec, "conduit", path_replacements)
        cfg.write("# conduit from uberenv\n")
        cfg.write(cmake_cache_entry("CONDUIT_DIR",conduit_dir))

        mfem_dir = get_spec_path(spec, "mfem", path_replacements)
        cfg.write("# mfem from uberenv\n")
        cfg.write(cmake_cache_entry("MFEM_DIR",mfem_dir))

        # optional tpls

        if "python" in spec or "devtools" in spec:
            python_bin_dir = get_spec_path(spec, "python", path_replacements, use_bin=True)
            cfg.write("# python from uberenv\n")
            cfg.write(cmake_cache_entry("PYTHON_EXECUTABLE",pjoin(python_bin_dir, "python")))
        else:
            cfg.write("# python not built by uberenv\n\n")

        # optional tpls (dev tools)

        if "doxygen" in spec:
            doxygen_bin_dir = get_spec_path(spec, "doxygen", path_replacements, use_bin=True)
            cfg.write("# doxygen from uberenv\n")
            cfg.write(cmake_cache_entry("DOXYGEN_EXECUTABLE", pjoin(doxygen_bin_dir, "doxygen")))
        else:
            cfg.write("# doxygen not built by uberenv\n\n")

        if "py-sphinx" in spec:
            python_bin_dir = get_spec_path(spec, "python", path_replacements, use_bin=True)
            cfg.write("# sphinx {} from uberenv\n".format(spec["py-sphinx"].version))
            cfg.write(cmake_cache_entry("SPHINX_EXECUTABLE", pjoin(python_bin_dir, "sphinx-build")))
        else:
            cfg.write("# sphinx not built by uberenv\n\n")

        if "py-shroud" in spec:
            python_bin_dir = get_spec_path(spec, "python", path_replacements, use_bin=True)
            cfg.write("# shroud {} from uberenv\n".format(spec["py-shroud"].version))
            cfg.write(cmake_cache_entry("SHROUD_EXECUTABLE", pjoin(python_bin_dir, "shroud")))
        else:
            cfg.write("# shroud not built by uberenv\n\n")

        if "uncrustify" in spec:
            uncrustify_bin_dir = get_spec_path(spec, "uncrustify", path_replacements, use_bin=True)
            cfg.write("# uncrustify from uberenv\n")
            cfg.write(cmake_cache_entry("UNCRUSTIFY_EXECUTABLE", pjoin(uncrustify_bin_dir, "uncrustify")))
        else:
            cfg.write("# uncrustify not built by uberenv\n\n")

        if "lcov" in spec:
            lcov_dir = get_spec_path(spec, "lcov", path_replacements)
            cfg.write("# lcov and genhtml from uberenv\n")
            cfg.write(cmake_cache_entry("LCOV_PATH", pjoin(lcov_dir,"usr","bin","lcov")))
            cfg.write(cmake_cache_entry("GENHTML_PATH",pjoin(lcov_dir,"usr","bin","genhtml")))
        else:
            cfg.write("# lcov and genhtml not built by uberenv\n\n")

        #MPI SUPPORT (when mpich is enabled)
        if "mpich" in spec:
            mpiexec  = pjoin(spec['mpich'].prefix.bin,"mpiexec")
            mpicc    = pjoin(spec['mpich'].prefix.bin,"mpicc")
            mpif90   = pjoin(spec['mpich'].prefix.bin,"mpif90")
            cfg.write("# MPI Support\n")
            cfg.write(cmake_cache_option("ENABLE_MPI",True))
            cfg.write(cmake_cache_entry("MPI_C_COMPILER",mpicc))
            # we use `mpicc` as `MPI_CXX_COMPILER` b/c we don't want to introduce 
            # linking deps to the MPI C++ libs (we aren't using C++ features of MPI)
            cfg.write(cmake_cache_entry("MPI_CXX_COMPILER",mpicc))
            cfg.write(cmake_cache_entry("MPI_Fortran_COMPILER",mpif90))
            cfg.write(cmake_cache_entry("MPIEXEC",mpiexec))

        cfg.write("##################################\n")
        cfg.write("# end uberenv host-config\n")
        cfg.write("##################################\n")
        cfg.write("\n")
        cfg.close()
        mkdirp(prefix)
        install(host_cfg_fname,prefix)
        print "[result host-config file: %s]" % host_cfg_fname

