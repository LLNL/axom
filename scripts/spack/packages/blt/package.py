# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.pkg.llnl.radiuss.blt import Blt as RadiussBlt

class Blt(RadiussBlt):
    # Note: This just points at a commit in the bugfix/white238/openmp branch.
    # It also has a made up version
    version("0.6.1.3", commit="1b6ababd46d3b1c1fcfa8b713c76121800be9c3e")
