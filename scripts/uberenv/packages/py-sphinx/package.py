from spack import *

class PySphinx(Package):
    """Sphinx Documentation Generator."""

    homepage = "http://sphinx-doc.org/"
    url      = "https://pypi.python.org/packages/source/S/Sphinx/Sphinx-1.3.1.tar.gz#md5=8786a194acf9673464c5455b11fd4332"

    version('1.3.1', '8786a194acf9673464c5455b11fd4332')

    depends_on("py-setuptools")

    extends('python')

    def install(self, spec, prefix):
        # sphinx doesn't play well with --prefix installs, for now simply 
        # install to the spack python
        python('setup.py', 'install') #, '--prefix=%s' % prefix)
