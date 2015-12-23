from spack import *

class PyPybindgen(Package):
    """Pybindgen is a tool that generates Python bindings for C/C++ code for Python."""

    homepage = "https://launchpad.net/pybindgen"
    url      = "https://pypi.python.org/packages/source/P/PyBindGen/PyBindGen-0.17.0.tar.gz"

    version('0.17.0', '7d8fe2b3b4646c3c1d9e5342b1645f6a')

    extends('python')

    def install(self, spec, prefix):
        # simply install to the spack python
        python('setup.py', 'install') 
