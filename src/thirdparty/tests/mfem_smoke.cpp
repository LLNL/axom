/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

//-----------------------------------------------------------------------------
///
/// file: mfem_smoke.cpp
///
//-----------------------------------------------------------------------------

#include "mfem.hpp"
#include "gtest/gtest.h"


//-----------------------------------------------------------------------------
TEST(mfem_smoke, basic_use)
{
  // Simple usage of a basic mfem type
  mfem::Element* el = new mfem::Quadrilateral(0,1,2,3);

  EXPECT_EQ( mfem::Element::QUADRILATERAL, el->GetType() );
  EXPECT_EQ( 4, el->GetNVertices() );

  delete el;
}
