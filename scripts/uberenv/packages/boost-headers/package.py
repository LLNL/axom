from spack import *
import shutil

import os
from os.path import join as pjoin

class BoostHeaders(Package):
    """Boost provides free peer-reviewed portable C++ source
       libraries, emphasizing libraries that work well with the C++
       Standard Library.

       Boost libraries are intended to be widely useful, and usable
       across a broad spectrum of applications. The Boost license
       encourages both commercial and non-commercial use.
    """
    homepage = "http://www.boost.org"
    url      = "http://downloads.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.bz2"
    list_url = "http://sourceforge.net/projects/boost/files/boost/"
    list_depth = 2

    version('1.58.0', 'b8839650e61e9c1c0a89f371dd475546')    

    def url_for_version(self, version):
        """Handle Boost's weird URLs, which write the version two different ways."""
        parts = [str(p) for p in Version(version)]
        dots = ".".join(parts)
        underscores = "_".join(parts)
        return "http://downloads.sourceforge.net/project/boost/boost/%s/boost_%s.tar.bz2" % (
            dots, underscores)

    def install(self, spec, prefix):
        # simply install the headers
        mkdirp(prefix.include)
        shutil.copytree('boost', pjoin(prefix.include,"boost"))
