from spack import *

class PyBreathe(Package):
    """Breathe provides a bridge between the Sphinx and Doxygen documentation systems."""

    homepage = "https://breathe.readthedocs.org/en/latest/"
    url      = "https://pypi.python.org/packages/source/b/breathe/breathe-4.0.0.tar.gz"

    version('4.0.0', '32316d5a890a3124ea3e8a9e0b2b3b97')

    extends('python')


    def install(self, spec, prefix):
        # breathe + sphinx doesn't play well with --prefix installs, for now simply 
        # install to the spack python
        python('setup.py', 'install') #, '--prefix=%s' % prefix)
