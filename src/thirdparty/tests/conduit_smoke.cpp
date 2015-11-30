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
/// file: conduit_smoke.cpp
///
//-----------------------------------------------------------------------------

#include "conduit.hpp"

#include <iostream>
#include "gtest/gtest.h"

//-----------------------------------------------------------------------------
TEST(conduit_smoke, basic_use)
{
    EXPECT_EQ(sizeof(conduit::uint32),(size_t)4);
    EXPECT_EQ(sizeof(conduit::uint64),(size_t)8);
    EXPECT_EQ(sizeof(conduit::float64),(size_t)8);
    
    std::cout << conduit::about() << std::endl;
    
}
