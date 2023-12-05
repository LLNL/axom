// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: mfem_smoke.cpp
///
//-----------------------------------------------------------------------------

#include "mfem.hpp"
#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(mfem_smoke, check_version)
{
  std::cout << "Using mfem version: "     //
            << MFEM_VERSION_MAJOR << "."  //
            << MFEM_VERSION_MINOR << "."  //
            << MFEM_VERSION_PATCH << std::endl;

  EXPECT_TRUE(MFEM_VERSION_MAJOR >= 0);
  EXPECT_TRUE(MFEM_VERSION_MINOR >= 0);
  EXPECT_TRUE(MFEM_VERSION_PATCH >= 0);
}

//-----------------------------------------------------------------------------
TEST(mfem_smoke, basic_use)
{
  // Simple usage of a basic mfem type
  mfem::Element* el = new mfem::Quadrilateral(0, 1, 2, 3);

  EXPECT_EQ(mfem::Element::QUADRILATERAL, el->GetType());
  EXPECT_EQ(4, el->GetNVertices());

  delete el;
}
