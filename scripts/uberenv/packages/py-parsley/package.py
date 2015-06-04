from spack import *

class PyParsley(Package):
    """Parsley is a parsing library for people who find parsers scary or annoying."""

    homepage = "https://pypi.python.org/pypi/Parsley"
    url      = "https://pypi.python.org/packages/source/P/Parsley/Parsley-1.2.tar.gz"

    version('1.2', '9905aa78f3256604f56e7bdd6aaf6123')

    extends('python')

    def install(self, spec, prefix):
        # simply install to the spack python
        python('setup.py', 'install') 