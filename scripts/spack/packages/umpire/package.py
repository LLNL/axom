# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.pkg.llnl.radiuss.umpire import Umpire as RadiussUmpire

class Umpire(RadiussUmpire):
    # Note: This just points at a commit in the task/update-blt-tpl-exports branch.
    # It also has a made up version
    version("2024.01.31", commit="f2ae5edd1e49e7fa18656203d742de3bd9490961", submodules=False)
