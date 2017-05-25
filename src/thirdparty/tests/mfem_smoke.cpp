/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
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
