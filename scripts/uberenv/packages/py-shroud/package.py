from spack import *
from spack.hooks.sbang import filter_shebang

class PyShroud(Package):
    """Create Fortran wrappers for a C++ library."""

    homepage = "https://github.com/LLNL/shroud"
    url      = "https://github.com/LLNL/shroud/archive/v0.10.0.tar.gz"

    version('0.10.0', 'c76bc9b1228d53ed0aabeb5806a486c7')
    version('0.9.0', 'f6a5ce0602a0c2d1d47de78c04ab302c')
    version('0.8.0', 'ec94d6f9cf3246d4370007abd4d270d8')

    extends('python')

    depends_on("py-pyyaml")

    def install(self, spec, prefix):
        # simply install to the spack python
        python('setup.py', 'install') 

        # shroud lives in python's bin dir
        shroud_scripts = ["shroud"]
        for script in shroud_scripts:
            script_path = join_path(spec["python"].prefix,"bin",script)
            # use spack sbang to fix issues with shebang that is too long
            filter_shebang(script_path)
