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

// This value is such that the 64Kb limit on device constant memory is not hit
// in check_alloc_realloc_free when reallocating to 3 * SIZE.
constexpr int SIZE = 5345;

class CopyTest :
  public ::testing::TestWithParam<::testing::tuple<std::string, std::string> >
{
public:
  void SetUp() override
  {
#ifdef AXOM_USE_UMPIRE
    umpire::ResourceManager & rm = umpire::ResourceManager::getInstance();
#endif

    src_string = ::testing::get<0>(GetParam());
    dst_string = ::testing::get<1>(GetParam());

    if (src_string == "NEW")
    {
      src_array = new int[SIZE];
    }
    else if (src_string == "MALLOC")
    {
      src_array = static_cast<int*>(std::malloc(SIZE * sizeof(int)));
    }
    else if (src_string == "STATIC")
    {
      src_array = m_static_src_array;
    }
#ifdef AXOM_USE_UMPIRE
    else
    {
      auto source_allocator = rm.getAllocator(src_string);
      src_array =
        static_cast<int*>(source_allocator.allocate(SIZE * sizeof(int)));
    }
#endif

    if (dst_string == "NEW")
    {
      dst_array = new int[SIZE];
    }
    else if (dst_string == "MALLOC")
    {
      dst_array = static_cast<int*>(std::malloc(SIZE * sizeof(int)));
    }
    else if (dst_string == "STATIC")
    {
      dst_array = m_static_dst_array;
    }
#ifdef AXOM_USE_UMPIRE
    else
    {
      auto source_allocator = rm.getAllocator(dst_string);
      dst_array =
        static_cast<int*>(source_allocator.allocate(SIZE * sizeof(int)));
    }
#endif
  }

  void TearDown() override
  {
#ifdef AXOM_USE_UMPIRE
    auto& rm = umpire::ResourceManager::getInstance();
#endif

    if (src_string == "NEW")
    {
      delete[] src_array;
    }
    else if (src_string == "MALLOC")
    {
      std::free(src_array);
    }
#ifdef AXOM_USE_UMPIRE
    else if (src_string != "STATIC")
    {
      rm.deallocate(src_array);
    }
#endif

    if (dst_string == "NEW")
    {
      delete[] dst_array;
    }
    else if (dst_string == "MALLOC")
    {
      std::free(dst_array);
    }
#ifdef AXOM_USE_UMPIRE
    else if (dst_string != "STATIC")
    {
      rm.deallocate(dst_array);
    }
#endif
  }

  std::string src_string;
  std::string dst_string;

  int* src_array = nullptr;
  int* dst_array = nullptr;
  int host_array[SIZE];

private:
  int m_static_src_array[SIZE];
  int m_static_dst_array[SIZE];
};

#ifdef AXOM_USE_UMPIRE
void check_alloc_and_free(
  umpire::Allocator allocator=
    axom::getAllocator(axom::getResourceAllocatorID(umpire::resource::Host)),
  bool hostAccessible=true  )
#else
void check_alloc_and_free( bool hostAccessible=true)
#endif
{
  for ( int size = 0 ; size <= SIZE ; size = size * 2 + 1 )
  {
#ifdef AXOM_USE_UMPIRE
    int* buffer = axom::allocate< int >( size, allocator );

    if (size > 0)
    {
      umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
      EXPECT_EQ( allocator.getId(), rm.getAllocator( buffer ).getId() );
    }
#else
    int* buffer = axom::allocate< int >( size );
#endif

    if ( hostAccessible )
    {
      for ( int i = 0 ; i < size ; ++i )
      {
        buffer[ i ] = i;
      }

      for ( int i = 0 ; i < size ; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    axom::deallocate( buffer );
    EXPECT_TRUE( buffer == nullptr );
  }
}

#ifdef AXOM_USE_UMPIRE
void check_alloc_realloc_free(
  umpire::Allocator allocator=
    axom::getAllocator(axom::getResourceAllocatorID(umpire::resource::Host)),
  bool hostAccessible=true )
#else
void check_alloc_realloc_free( bool hostAccessible=true )
#endif
{
  for ( int size = 0 ; size <= SIZE ; size = size * 2 + 1 )
  {
    int buffer_size = size;

#ifdef AXOM_USE_UMPIRE
    int* buffer = axom::allocate< int >( buffer_size, allocator );

    umpire::ResourceManager & rm = umpire::ResourceManager::getInstance();
    if (buffer_size > 0)
    {
      ASSERT_EQ(allocator.getId(), rm.getAllocator(buffer).getId());
    }
#else
    int* buffer = axom::allocate< int >( buffer_size );
#endif

    if ( hostAccessible )
    {
      // Populate the buffer.
      for ( int i = 0 ; i < buffer_size ; ++i )
      {
        buffer[ i ] = i;
      }

      // Check the values.
      for ( int i = 0 ; i < buffer_size ; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    // Reallocate to a larger size.
    buffer_size *= 3;
    buffer = axom::reallocate< int >( buffer, buffer_size );
#ifdef AXOM_USE_UMPIRE
    if (buffer_size > 0)
    {
      ASSERT_EQ(allocator.getId(), rm.getAllocator(buffer).getId());
    }
#endif

    if ( hostAccessible )
    {
      // Populate the new values.
      for ( int i = size ; i < buffer_size ; ++i )
      {
        buffer[ i ] = i;
      }

      // Check all the values.
      for ( int i = 0 ; i < buffer_size ; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    // Reallocate to a smaller size.
    buffer_size /= 5;
    buffer = axom::reallocate< int >( buffer, buffer_size );
#ifdef AXOM_USE_UMPIRE
    if (buffer_size > 0)
    {
      ASSERT_EQ(allocator.getId(), rm.getAllocator(buffer).getId());
    }
#endif

    if ( hostAccessible )
    {
      // Check all the values.
      for ( int i = 0 ; i < buffer_size ; ++i )
      {
        EXPECT_EQ( buffer[ i ], i );
      }
    }

    // Free
    axom::deallocate( buffer );
    EXPECT_TRUE( buffer == nullptr );
  }
}

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

#ifdef AXOM_USE_UMPIRE

TEST( core_memory_management, set_get_default_memory_space )
{
  const int HostAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Host );
  EXPECT_EQ( HostAllocatorID, axom::getDefaultAllocator().getId() );

#ifdef AXOM_USE_CUDA
  const int PinnedAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Pinned );

  axom::setDefaultAllocator( axom::getAllocator( PinnedAllocatorID ) );
  EXPECT_EQ( PinnedAllocatorID, axom::getDefaultAllocator().getId() );

  const int DeviceAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Device );
  axom::setDefaultAllocator( axom::getAllocator( DeviceAllocatorID ) );
  EXPECT_EQ( DeviceAllocatorID, axom::getDefaultAllocator().getId() );

  const int ConstantAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Constant );
  axom::setDefaultAllocator( axom::getAllocator( ConstantAllocatorID ) );
  EXPECT_EQ( ConstantAllocatorID, axom::getDefaultAllocator().getId() );

  const int UnifiedAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Unified );
  axom::setDefaultAllocator( axom::getAllocator( UnifiedAllocatorID ) );
  EXPECT_EQ( UnifiedAllocatorID, axom::getDefaultAllocator().getId() );
#endif


  axom::setDefaultAllocator( axom::getAllocator( HostAllocatorID ) );
  EXPECT_EQ( HostAllocatorID, axom::getDefaultAllocator().getId() );
}

#endif

//------------------------------------------------------------------------------
TEST( core_memory_management, alloc_free )
{
#ifdef AXOM_USE_UMPIRE

  constexpr bool HOST_ACCESSIBLE = true;

  const int HostAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Host );
  check_alloc_and_free( axom::getAllocator( HostAllocatorID ),
                        HOST_ACCESSIBLE
                        );

#ifdef AXOM_USE_CUDA

  constexpr bool NOT_HOST_ACCESSIBLE = false;

  const int PinnedAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Pinned );
  check_alloc_and_free( axom::getAllocator( PinnedAllocatorID ),
                        HOST_ACCESSIBLE
                        );

  const int DeviceAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Device );
  check_alloc_and_free( axom::getAllocator( DeviceAllocatorID ),
                        NOT_HOST_ACCESSIBLE
                        );

  const int ConstantAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Constant );
  check_alloc_and_free( axom::getAllocator( ConstantAllocatorID ),
                        NOT_HOST_ACCESSIBLE
                        );

  const int UnifiedAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Unified );
  check_alloc_and_free( axom::getAllocator( UnifiedAllocatorID ),
                        HOST_ACCESSIBLE
                        );
#endif

#endif

  check_alloc_and_free();
}

//------------------------------------------------------------------------------
TEST( core_memory_management, alloc_realloc_free )
{
#ifdef AXOM_USE_UMPIRE

  constexpr bool HOST_ACCESSIBLE = true;

  check_alloc_realloc_free( axom::getAllocator( umpire::resource::Host ),
                            HOST_ACCESSIBLE
                            );

#ifdef AXOM_USE_CUDA

  constexpr bool NOT_HOST_ACCESSIBLE = false;

  const int PinnedAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Pinned );
  check_alloc_realloc_free( axom::getAllocator( PinnedAllocatorID ),
                            HOST_ACCESSIBLE
                            );

  const int DeviceAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Device );
  check_alloc_realloc_free( axom::getAllocator( DeviceAllocatorID ),
                            NOT_HOST_ACCESSIBLE
                            );

  // Umpire doesn't allow reallocation of Constant memory.
  // check_alloc_realloc_free( axom::getAllocator( umpire::resource::Constant ),
  // false );

  const int UnifiedAllocatorID =
    axom::getResourceAllocatorID( umpire::resource::Unified );
  check_alloc_realloc_free( axom::getAllocator( UnifiedAllocatorID ),
                            HOST_ACCESSIBLE
                            );
#endif

#endif

  check_alloc_realloc_free();
}

TEST_P(CopyTest, Copy)
{
  std::cout << "SRC = " << src_string << ", DST = " << dst_string << std::endl;
  for (int i = 0 ; i < SIZE ; ++i)
  {
    host_array[i] = i;
  }

  axom::copy(src_array, host_array, SIZE * sizeof(int));
  axom::copy(dst_array, src_array, SIZE * sizeof(int));

  for (int i = 0 ; i < SIZE ; ++i)
  {
    host_array[i] = -i;
  }

  axom::copy(host_array, src_array, SIZE * sizeof(int));

  for (int i = 0 ; i < SIZE ; ++i)
  {
    ASSERT_EQ(host_array[i], i);
  }
}

const std::string copy_locations[] = {
  "NEW", "MALLOC", "STATIC"
#if defined(AXOM_USE_UMPIRE)
  , "HOST"
#if defined(UMPIRE_ENABLE_DEVICE)
  , "DEVICE"
#endif
#if defined(UMPIRE_ENABLE_UM)
  , "UM"
#endif
#if defined(UMPIRE_ENABLE_PINNED)
  , "PINNED"
#endif
#endif
};


INSTANTIATE_TEST_CASE_P(
  core_memory_management,
  CopyTest,
  ::testing::Combine(
    ::testing::ValuesIn(copy_locations),
    ::testing::ValuesIn(copy_locations)
    ));

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int result = RUN_ALL_TESTS();
  return( result );
}
