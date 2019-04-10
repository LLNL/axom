// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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

constexpr int ALLOCATION_SIZE = 5345;

#ifdef AXOM_USE_UMPIRE
void check_alloc_and_free( umpire::Allocator allocator=axom::getAllocator( umpire::resource::Host ), bool hostAccessable=true  )
#else
void check_alloc_and_free( bool hostAccessable=true)
#endif
{
  for ( int size = 0; size <= ALLOCATION_SIZE; size = size * 2 + 1 )
  {
#ifdef AXOM_USE_UMPIRE
    int* buffer = axom::alloc< int >( size, allocator );
    
    if (size > 0)
    {
      umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
      umpire::Allocator allocator_cpy = rm.getAllocator( buffer );
      EXPECT_EQ( allocator_cpy.getId(), allocator.getId() );
    }
#else
    int* buffer = axom::alloc< int >( size );
#endif

    if ( hostAccessable )
    {
      for ( int i = 0; i < size; ++i )
      {
        buffer[ i ] = i;
      }

      for ( int i = 0; i < size; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    axom::free( buffer );
    EXPECT_TRUE( buffer == nullptr );
  }
}

#ifdef AXOM_USE_UMPIRE
void check_alloc_realloc_free( umpire::Allocator allocator=axom::getAllocator( umpire::resource::Host ), bool hostAccessable=true )
#else
void check_alloc_realloc_free( bool hostAccessable=true )
#endif
{
  for ( int size = 0; size <= ALLOCATION_SIZE; size = size * 2 + 1 )
  {
    int buffer_size = size;

#ifdef AXOM_USE_UMPIRE
    int* buffer = axom::alloc< int >( buffer_size, allocator );
#else
    int* buffer = axom::alloc< int >( buffer_size );
#endif

    if ( hostAccessable )
    {
      // Populate the buffer.
      for ( int i = 0; i < buffer_size; ++i )
      {
        buffer[ i ] = i;
      }

      // Check the values.
      for ( int i = 0; i < buffer_size; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    // Reallocate to a larger size.
    buffer_size *= 3;
    buffer = axom::realloc< int >( buffer, buffer_size );

    if ( hostAccessable )
    {
      // Populate the new values.
      for ( int i = size; i < buffer_size; ++i )
      {
        buffer[ i ] = i;
      }

      // Check all the values.
      for ( int i = 0; i < buffer_size; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    // Reallocate to a smaller size.
    buffer_size /= 5;
    buffer = axom::realloc< int >( buffer, buffer_size );

    if ( hostAccessable )
    {
      // Check all the values.
      for ( int i = 0; i < buffer_size; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    // Free
    axom::free( buffer );
    EXPECT_TRUE( buffer == nullptr );
  }
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

#ifdef AXOM_USE_UMPIRE

TEST( core_memory_management, set_get_default_memory_space )
{
  EXPECT_EQ( umpire::resource::Host, axom::getDefaultAllocator().getId() );

#ifdef AXOM_USE_CUDA
  axom::setDefaultAllocator( axom::getAllocator( umpire::resource::Pinned ) );
  EXPECT_EQ( umpire::resource::Pinned, axom::getDefaultAllocator().getId() );

  axom::setDefaultAllocator( axom::getAllocator( umpire::resource::Device ) );
  EXPECT_EQ( umpire::resource::Device, axom::getDefaultAllocator().getId() );

  axom::setDefaultAllocator( axom::getAllocator( umpire::resource::Constant ) );
  EXPECT_EQ( umpire::resource::Constant, axom::getDefaultAllocator().getId() );

  axom::setDefaultAllocator( axom::getAllocator( umpire::resource::Unified ) );
  EXPECT_EQ( umpire::resource::Unified, axom::getDefaultAllocator().getId() );
#endif

  axom::setDefaultAllocator( axom::getAllocator( umpire::resource::Host ) );
}

#endif

//------------------------------------------------------------------------------
TEST( core_memory_management, alloc_free )
{
#ifdef AXOM_USE_UMPIRE

  check_alloc_and_free( axom::getAllocator( umpire::resource::Host ), true );

#ifdef AXOM_USE_CUDA
  check_alloc_and_free( axom::getAllocator( umpire::resource::Pinned ), true );
  check_alloc_and_free( axom::getAllocator( umpire::resource::Device ), false );
  check_alloc_and_free( axom::getAllocator( umpire::resource::Constant ), false );
  check_alloc_and_free( axom::getAllocator( umpire::resource::Unified ), true );
#endif

#endif

  check_alloc_and_free();
}

//------------------------------------------------------------------------------
TEST( core_memory_management, alloc_realloc_free )
{
#ifdef AXOM_USE_UMPIRE

  check_alloc_realloc_free( axom::getAllocator( umpire::resource::Host ), true );

#ifdef AXOM_USE_CUDA
  check_alloc_realloc_free( axom::getAllocator( umpire::resource::Pinned ), true );
  check_alloc_realloc_free( axom::getAllocator( umpire::resource::Device ), false );
  // Umpire doesn't allow reallocation of Constant memory.
  // check_alloc_realloc_free( axom::getAllocator( umpire::resource::Constant ), false );
  check_alloc_realloc_free( axom::getAllocator( umpire::resource::Unified ), true );
#endif

#endif

  check_alloc_realloc_free();
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int result = RUN_ALL_TESTS();
  return( result );
}
