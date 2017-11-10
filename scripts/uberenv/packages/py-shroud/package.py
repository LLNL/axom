from spack import *
from spack.hooks.sbang import filter_shebang

class PyShroud(Package):
    """Create Fortran wrappers for a C++ library."""

    homepage = "https://github.com/LLNL/shroud"
    url      = "https://github.com/LLNL/shroud/archive/v0.2.tar.gz"

    version('0.2', 'af7b53240c20dcdace0ee0ea8e3e9818')

    extends('python')

    depends_on("py-pyyaml")
    depends_on("py-parsley")

    def install(self, spec, prefix):
        # simply install to the spack python
        python('setup.py', 'install') 

        # shroud lives in python's bin dir
        shroud_scripts = ["shroud"]
        for script in shroud_scripts:
            script_path = join_path(spec["python"].prefix,"bin",script)
            # use spack sbang to fix issues with shebang that is too long
            filter_shebang(script_path)
