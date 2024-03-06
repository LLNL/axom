# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.pkg.llnl.radiuss.blt import Blt as RadiussBlt

class Blt(RadiussBlt):
    # Note: This just points at a commit in the bugfix/white238/openmp branch.
    # It also has a made up version
    version("0.6.1.4", commit="c98f320835d71f778fbffcc48f07675142c08635")
