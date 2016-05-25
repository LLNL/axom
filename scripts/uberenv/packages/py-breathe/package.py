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
        # sphinx_build lives in python's bin dir and for some reason breathe's
        # setuptools install reinstalls the sphinx build scripts (even though
        # it finds the existing sphinx install).
        # To keep this from undermining us, we need to patch the new copies
        # of these scripts.
        sphinx_scripts = ["sphinx-apidoc",
                          "sphinx-autogen",
                          "sphinx-build",
                          "sphinx-quickstart"]
        for script in sphinx_scripts:
            script_path = join_path(spec["python"].prefix,"bin",script)
            # use spack sbang to fix issues with shebang that is too long
            filter_shebang(script_path)