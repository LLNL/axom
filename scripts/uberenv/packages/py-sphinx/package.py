from spack import *
from spack.hooks.sbang import filter_shebang


class PySphinx(Package):
    """Sphinx Documentation Generator."""

    homepage = "http://sphinx-doc.org/"
    url      = "https://pypi.python.org/packages/source/S/Sphinx/Sphinx-1.3.6.tar.gz#md5=7df638f47749f9284889c93012ffa07f"

    version('1.3.6', '7df638f47749f9284889c93012ffa07f')

    depends_on("py-setuptools")
    depends_on("py-six")
    depends_on("py-jinja2")
    depends_on("py-pygments")
    depends_on("py-docutils")
    depends_on("py-snowballstemmer")
    depends_on("py-babel")
    depends_on("py-alabaster")
    depends_on("py-sphinx-rtd-theme")

    extends('python')

    def install(self, spec, prefix):
        # sphinx doesn't play well with --prefix installs, for now simply 
        # install to the spack python
        python('setup.py', 'install') #, '--prefix=%s' % prefix)
        # sphinx_build lives in python's bin dir
        sphinx_scripts = ["sphinx-apidoc",
                          "sphinx-autogen",
                          "sphinx-build",
                          "sphinx-quickstart"]
        for script in sphinx_scripts:
            script_path = join_path(spec["python"].prefix,"bin",script)
            # use spack sbang to fix issues with shebang that is too long
            filter_shebang(script_path)
        

