##############################################################################
# Copyright (c) 2013, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This is a configuration file for building lcov using Spack.
# -- Aaron Black black27@llnl.gov 7/6/2015
#
# For details on Spack:
# See https://scalability-llnl.github.io/spack
#
##############################################################################
from spack import *

class Lcov(Package):
    """LCOV is a graphical front-end for GCC's coverage testing tool gcov. It
       collects gcov data for multiple source files and creates HTML pages
       containing the source code annotated with coverage information. It also
       adds overview pages for easy navigation within the file structure. LCOV
       supports statement, function and branch coverage measurement."""
    homepage  = 'https://ltp.sourceforget.net/coverage/lcov.php'

    version('1.11', 'e79b799ae3ce149aa924c7520e993024',
            url = 'http://downloads.sourceforge.net/ltp/lcov-1.11.tar.gz')

    def install(self, spec, prefix):
        make('install', 'PREFIX=%s' % prefix)
