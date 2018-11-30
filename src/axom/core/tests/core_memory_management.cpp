/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

#include "gtest/gtest.h"

#include "axom/core/memory_management.hpp"

TEST( core_memory_managment, allocation )
{
  std::cout<<"Testing allocation functions."<< std::endl;
  constexpr int MAX_SIZE = 1048576 ;
  for ( int initial_size = 2; initial_size <= MAX_SIZE; initial_size *= 2 )
  {
    int buffer_size = initial_size;
    int* buffer = axom::alloc< int >( buffer_size );

    for (int i = 0 ; i < buffer_size ; i++)
    {
      buffer[i] = i;
    }

    buffer_size *= 2;
    buffer = axom::realloc( buffer, buffer_size );
    for (int i = 0 ; i < buffer_size / 2 ; i++)
    {
      EXPECT_EQ( buffer[i], i );
    }

    buffer_size /= 4;
    buffer = axom::realloc( buffer, buffer_size );
    for (int i = 0 ; i < buffer_size ; i++)
    {
      EXPECT_EQ( buffer[i], i );
    }

    axom::free( buffer );
    EXPECT_TRUE( buffer == nullptr );
  }
}
