###############################################################################
# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
###############################################################################
from spack import *

class PyBreathe(PythonPackage):
    """Breathe provides a bridge between the Sphinx and Doxygen documentation systems."""

    homepage = "https://breathe.readthedocs.org/en/latest/"
    url      = "https://pypi.python.org/packages/source/b/breathe/breathe-4.0.0.tar.gz"

    version('4.0.0', '32316d5a890a3124ea3e8a9e0b2b3b97')

    extends('python', ignore='bin/(pybabel|pygmentize)')

    # See here for upstream list of dependencies:
    # https://github.com/michaeljones/breathe/blob/master/setup.py
    
    depends_on("python@2.7:2.8,3.4", type=('build','run'))
    depends_on("py-setuptools",      type=('build','run'))
    depends_on("py-sphinx@1.4:",     type=('build','run'))
    depends_on('py-docutils@0.5:',   type=('build','run'))
    depends_on('py-six@1.4:',        type=('build','run'))

