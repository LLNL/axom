# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class Axomdevtools(BundlePackage):
    """This is a set of tools necessary for the developers of Axom"""

    version('fakeversion')

    maintainers = ['white238']

    depends_on("python")
    depends_on("doxygen")
    depends_on("cppcheck")
    depends_on("graphviz")
    depends_on("py-sphinx")
    depends_on("py-shroud")
    depends_on("llvm+clang@10.0.0")
