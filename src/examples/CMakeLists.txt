# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

###############################################################################
#
# file: src/examples/CMakeLists.txt
#
###############################################################################

# install using examples
install(DIRECTORY using-with-cmake using-with-make
        DESTINATION examples/axom)
# note: this is installed to examples/axom to avoid name conflicts with other
# libs (this will crop up if we use spack to install axom as a python 
# extension)
    

