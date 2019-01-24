from spack import *

import socket
import os
import platform
from os.path import join as pjoin

def cmake_cache_entry(name, value, comment=""):
    """Generate a string for a cmake cache variable""" 
    return 'set(%s "%s" CACHE PATH "%s")\n\n' % (name,value,comment)

def cmake_cache_option(name, boolean_value, comment=""):
    """Generate a string for a cmake configuration option""" 
    
    value = "ON" if boolean_value else "OFF"
    return 'set(%s %s CACHE BOOL "%s")\n\n' % (name,value,comment)

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
    variant('devtools', default=False, description="Build development tools (such as sphinx, uncrustify, etc)")

    # use ~cmake to skip cmake build and use whatever cmake is in the
    # users path 
    # (given the pain of building cmake on BGQ, this is really only for BGQ)
    variant('cmake',   default=True, description="Build cmake.")

    variant("python",   default=False, description="Build python")
    variant("mpi",      default=True, description="Build MPI support")
    variant("mfem",     default=True, description="Build mfem")

    variant("hdf5",     default=True, description="Build hdf5")

    variant("scr",      default=False, description="Build SCR")

    depends_on("hdf5~cxx~shared~fortran", when="+hdf5")

    depends_on("conduit~shared+python",when="+python")
    depends_on("conduit~shared~python",when="~python")
    depends_on("conduit~shared+hdf5+python",when="+hdf5+python")
    depends_on("conduit~shared~hdf5+python",when="~hdf5+python")
    depends_on("conduit~shared+hdf5~python",when="+hdf5~python")
    depends_on("conduit~shared~hdf5~python",when="~hdf5~python")
    
    depends_on("scr", when="+scr")

    # builds serial version of mfem that does not depend on Sidre
    depends_on("mfem~hypre~metis~mpi~gzstream",   when="+mfem")

    # optional tpl builds

    # NOTE: Due to spack defaulting to the newest version given in a package
    # The version of CMake is controlled in the <sys_type>/packages.yaml file.
    # This is a problem when multiple packages depend on another package and
    # the required (or not required) versions don't match
    depends_on("cmake", when="+cmake")

    if "darwin" in platform.system().lower():
        depends_on("mpich@3.0.4")
        depends_on("openssl@1.0.2j")

    depends_on("python",    when="+devtools")
    depends_on("doxygen",   when="+devtools")
    depends_on("uncrustify@0.61",when="+devtools")
    depends_on("cppcheck",when="+devtools")

    depends_on("python",   when="+python")

    depends_on("py-sphinx", when="+devtools")
    depends_on("py-breathe", when="+devtools")
    depends_on("py-shroud", when="+devtools")

    depends_on("mpi",when="+mpi")

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

        # are we on a specific machine 
        on_bgq = 'bgq' in sys_type
        on_blueos = 'blueos' in sys_type
        on_toss =  'toss_3' in sys_type

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
        cfg.write("# !!!! This is a generated file, edit at own risk !!!!\n")
        cfg.write("##################################\n\n")
        cfg.write("##################################\n")
        cfg.write("#\n")
        cfg.write("# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.\n")
        cfg.write("#\n")
        cfg.write("# Produced at the Lawrence Livermore National Laboratory.\n")
        cfg.write("#\n")
        cfg.write("# LLNL-CODE-741217\n")
        cfg.write("#\n")
        cfg.write("# All rights reserved.\n")
        cfg.write("#\n")
        cfg.write("# This file is part of Axom.\n")
        cfg.write("#\n")
        cfg.write("# For details about use and distribution, please read axom/LICENSE.\n")
        cfg.write("#\n")
        cfg.write("##################################\n\n")
        cfg.write("##################################\n\n")
        cfg.write("# SYS_TYPE: %s\n" % (sys_type))
        cfg.write("# Compiler Spec: %s\n" % (spec.compiler))
        cfg.write("##################################\n\n")
        # show path to cmake for reference and to be used by config-build.py
        cfg.write("# CMake executable path: %s\n\n" % cmake_exe)

        # compiler settings
        cfg.write("##############\n")
        cfg.write("# Compilers\n")
        cfg.write("##############\n\n")

        if on_bgq:
            cfg.write("# Note: we build TPLs with the serial compiler then use MPI wrappers on bgq\n")
            cfg.write("# Serial compilers used by spack:\n")
            cfg.write("# C compiler: {0}\n".format(c_compiler))
            cfg.write("# C++ compiler: {0}\n\n".format(cpp_compiler))
            cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", spec['mpi'].prefix.bin.mpiclang))
            cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", spec['mpi'].prefix.bin.mpiclang + "++"))
        else:
            cfg.write("# C compiler used by spack\n")
            cfg.write(cmake_cache_entry("CMAKE_C_COMPILER",c_compiler))
            cfg.write("# C++ compiler used by spack\n")
            cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER",cpp_compiler))

        if f_compiler is None:
            cfg.write("# No fortran compiler\n\n")
            cfg.write(cmake_cache_option("ENABLE_FORTRAN",False))
        else:
            cfg.write("# Fortran compiler used by spack\n")
            cfg.write(cmake_cache_option("ENABLE_FORTRAN",True))
            cfg.write(cmake_cache_entry("CMAKE_Fortran_COMPILER",f_compiler))

        # TPL locations
        cfg.write("##############\n")
        cfg.write("# TPLs\n")
        cfg.write("##############\n\n")

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

        if "doxygen" in spec or "py-sphinx" in spec:
            cfg.write(cmake_cache_option("ENABLE_DOCS", True))

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
        else:
            cfg.write(cmake_cache_option("ENABLE_DOCS", False))

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

        if "cppcheck" in spec:
            cppcheck_bin_dir = get_spec_path(spec, "cppcheck", path_replacements, use_bin=True)
            cfg.write("# cppcheck from uberenv\n")
            cfg.write(cmake_cache_entry("CPPCHECK_EXECUTABLE", pjoin(cppcheck_bin_dir, "cppcheck")))
        else:
            cfg.write("# cppcheck not built by uberenv\n\n")

        cfg.write("##############\n")
        cfg.write("# MPI\n")
        cfg.write("##############\n\n")

        if "+mpi" in spec:
            cfg.write(cmake_cache_option("ENABLE_MPI", True))
            if on_bgq:
                cfg.write(cmake_cache_option("ENABLE_FIND_MPI", False,
                   "Use wrapper directly to stop FindMPI returning the wrong linker flags."))
                cfg.write(cmake_cache_option("ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC", True,
                    "Ensures that tests will be wrapped with srun to run on the backend nodes"))
                cfg.write(cmake_cache_entry("MPI_Fortran_INCLUDE_PATH",
                    "/usr/local/tools/deg/drivers/V1R2M0/ppc64/comm/gcc/include",
                    "Pass in an explicit path to help find mpif.h"))
            else:
                cfg.write(cmake_cache_entry("MPI_C_COMPILER", spec['mpi'].mpicc))
                cfg.write(cmake_cache_entry("MPI_CXX_COMPILER",
                                            spec['mpi'].mpicxx))
                cfg.write(cmake_cache_entry("MPI_Fortran_COMPILER",
                                            spec['mpi'].mpifc))

            # Determine MPIEXEC
            if on_blueos:
                mpiexec = join_path(spec['mpi'].prefix.bin, 'mpirun')
            else:
                mpiexec = join_path(spec['mpi'].prefix.bin, 'mpiexec')
                if not os.path.isfile(mpiexec):
                    mpiexec = "/usr/bin/srun"
            # starting with cmake 3.10, FindMPI expects MPIEXEC_EXECUTABLE
            # vs the older versions which expect MPIEXEC
            if self.spec["cmake"].satisfies('@3.10:'):
                cfg.write(cmake_cache_entry("MPIEXEC_EXECUTABLE", mpiexec))
            else:
                cfg.write(cmake_cache_entry("MPIEXEC", mpiexec))

            # Determine MPIEXEC_NUMPROC_FLAG
            if on_blueos:
                cfg.write(cmake_cache_entry("MPIEXEC_NUMPROC_FLAG", "-np"))
                cfg.write(cmake_cache_entry("BLT_MPI_COMMAND_APPEND", "mpibind"))
            else:
                # TODO: see if spack has this
                cfg.write(cmake_cache_entry("MPIEXEC_NUMPROC_FLAG", "-n"))
        else:
            cfg.write(cmake_cache_option("ENABLE_MPI", False))

        ##################################
        # Other machine specifics
        ##################################

        cfg.write("##############\n")
        cfg.write("# Other machine specifics\n")
        cfg.write("##############\n\n")

        # Enable death tests everwhere but BGQ
        if on_bgq:
            cfg.write(cmake_cache_option("ENABLE_GTEST_DEATH_TESTS", False))
        else:
            cfg.write(cmake_cache_option("ENABLE_GTEST_DEATH_TESTS", True))

        # BGQ
        if on_bgq:
            if "xlf" in str(spec.compiler):
                cfg.write(cmake_cache_entry("BLT_FORTRAN_FLAGS", "-WF,-C!",
                    "Converts C-style comments to Fortran style in preprocessed files"))

            cfg.write("# Manually set up HDF5 library dependencies for BGQ to bypass errors from CMake's FindHDF5\n")
            cfg.write(cmake_cache_entry("HDF5_C_LIBRARY_m", "-lm"))
            cfg.write(cmake_cache_entry("HDF5_C_LIBRARY_dl", "-ldl"))

            cfg.write(cmake_cache_option("CMAKE_SKIP_RPATH", True))

        # BlueOS
        elif on_blueos:
            if "xlf" in f_compiler:
                cfg.write(cmake_cache_entry("CMAKE_Fortran_COMPILER_ID", "XL",
                    "All of BlueOS compilers report clang due to nvcc, override to proper compiler family"))
            if "xlc" in c_compiler:
                cfg.write(cmake_cache_entry("CMAKE_C_COMPILER_ID", "XL",
                    "All of BlueOS compilers report clang due to nvcc, override to proper compiler family"))
            if "xlC" in cpp_compiler:
                cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER_ID", "XL",
                    "All of BlueOS compilers report clang due to nvcc, override to proper compiler family"))

            if "clang@coral_xlf" == str(spec.compiler):
                cfg.write(cmake_cache_entry("BLT_FORTRAN_FLAGS", "-WF,-C!",
                    "Converts C-style comments to Fortran style in preprocessed files"))
                cfg.write(cmake_cache_entry("BLT_EXE_LINKER_FLAGS",
                    "-Wl,-rpath,/usr/tce/packages/xl/xl-2018.05.18/lib/", 
                    "Adds a missing rpath for libraries associated with the fortran compiler"))
            elif "xl@coral" == str(spec.compiler):
                cfg.write(cmake_cache_entry("BLT_FORTRAN_FLAGS", "-WF,-C! -qxlf2003=polymorphic",
                    "Convert C-style comments to Fortran and link fortran exes to C++ libraries"))

        # TOSS3
        elif on_toss:
            if "gcc@4.9.3" == str(spec.compiler):
                cfg.write(cmake_cache_entry("SCR_DIR",
                    "/usr/gapps/axom/thirdparty_libs/scr-1.2.1/toss_3_x86_64_ib/gcc-4.9.3"))

        cfg.write("\n")
        cfg.close()
        mkdirp(prefix)
        install(host_cfg_fname,prefix)
        print "[result host-config file: %s]" % host_cfg_fname

