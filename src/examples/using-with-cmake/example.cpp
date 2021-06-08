// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


//-----------------------------------------------------------------------------
///
/// file: example.cpp
///
//-----------------------------------------------------------------------------

#include "axom/core/utilities/About.hpp"
#include "fmt/fmt.hpp"

#include <iostream>

int main()
{
   // Using fmt library exported by axom
   std::cout << fmt::format(
        "Example of using and installed version of axom v{}.{}.{}-{}",
        AXOM_VERSION_MAJOR, AXOM_VERSION_MINOR,
        AXOM_VERSION_PATCH, AXOM_VERSION_EXTRA) << std::endl << std::endl;

   // Uses installed axom library
   axom::about();
}

