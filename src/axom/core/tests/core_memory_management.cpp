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

#ifdef AXOM_USE_UMPIRE
#include "umpire/config.hpp"
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#endif

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------

void check_alloc_and_free( axom::MemorySpace spaceId )
{
#ifdef AXOM_USE_UMPIRE

  constexpr int MAX_SIZE = 1048576;

  int* data = axom::alloc< int >( MAX_SIZE, spaceId );
  auto& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator allocator = rm.getAllocator( data );
  EXPECT_EQ( allocator.getId(), axom::internal::umpire_type[ spaceId ] );

  axom::free( data );
  EXPECT_TRUE( data == nullptr );

#else

  /* silence compiler warnings */
  static_cast< void >( spaceId );

#endif
}

//------------------------------------------------------------------------------
void check_alloc_realloc_free( axom::MemorySpace spaceId )
{
  std::cout<<"Testing allocation functions."<< std::endl;
  constexpr int MAX_SIZE = 1048576 ;

  for ( int initial_size = 2; initial_size <= MAX_SIZE; initial_size *= 2 )
  {
    int buffer_size = initial_size;
    int* buffer = axom::alloc< int >( buffer_size, spaceId );
    EXPECT_TRUE( buffer != nullptr );

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

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( core_memory_management, alloc_free )
{
  check_alloc_and_free( axom::HOST );

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  check_alloc_and_free( axom::HOST );
  check_alloc_and_free( axom::DEVICE );
  check_alloc_and_free( axom::DEVICE_CONSTANT );
  check_alloc_and_free( axom::UNIFIED_MEMORY );
#endif

}

//------------------------------------------------------------------------------
TEST( core_memory_management, alloc_realloc_free )
{
  check_alloc_realloc_free( axom::HOST );

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  check_alloc_realloc_free( axom::UNIFIED_MEMORY );
#endif
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int result = RUN_ALL_TESTS();
  return( result );
}
