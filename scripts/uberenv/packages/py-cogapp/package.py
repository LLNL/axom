from spack import *

class PyCogapp(Package):
    """A code generator for executing Python snippets in source files."""

    homepage = "http://nedbatchelder.com/code/cog"
    url      = "https://pypi.python.org/packages/source/c/cogapp/cogapp-2.4.tar.gz"

    version('2.4', 'ac927d494c3e7ab6a39e61b6253c69eb')

    extends('python')

    def install(self, spec, prefix):
        # simply install to the spack python
        python('setup.py', 'install') 
