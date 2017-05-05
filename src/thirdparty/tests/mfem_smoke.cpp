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
TEST(mfem_smoke, configuration)
{
    // Verify that this copy of mfem is built without MPI
    bool hasMPI = false;
    #ifdef MFEM_USE_MPI
       hasMPI = true;
    #endif
    EXPECT_FALSE(hasMPI) << "Axom expects mfem to be built without MPI";

    // Verify that this copy of mfem is built without Sidre
    bool hasSidre = false;
    #ifdef MFEM_USE_SIDRE
       hasSidre = true;
    #endif
    EXPECT_FALSE(hasSidre) << "Axom expects mfem to be built without Sidre";
}

//-----------------------------------------------------------------------------
TEST(mfem_smoke, basic_use)
{
    // Simple usage of a basic mfem type
    mfem::Element* el = new mfem::Quadrilateral(0,1,2,3);    
    
    EXPECT_EQ( mfem::Element::QUADRILATERAL, el->GetType() );
    EXPECT_EQ( 4, el->GetNVertices() );
    
    delete el;
}
