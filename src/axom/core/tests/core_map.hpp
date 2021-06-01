// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/StackArray.hpp" /* for axom::StackArray */
#include "axom/core/Map.hpp"
#include "gtest/gtest.h"            /* for TEST and EXPECT_* macros */
#include <string>

namespace axom
{

TEST(core_map, initialization)
{
  constexpr int N = 10; /* Number of values to store */
  Map <int, int> test(N);

}


} /* namespace axom */
