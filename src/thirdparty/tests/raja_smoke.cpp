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
#include "axom/config.hpp"                    // for compile-time definitions
#include "axom/core/utilities/Utilities.hpp"  // alloc() / free() methods

// RAJA includes
#include "RAJA/RAJA.hpp"     // for RAJA

// Google Test
#include "gtest/gtest.h"     // for google test functions

// C/C++ includes
#include <iostream> // for std::cout

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

//------------------------------------------------------------------------------
template < typename T >
T* allocate( std::size_t N )
{
  T* ptr = nullptr;
#if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA)
  ptr = cudaMallocManaged( (void**)&ptr, sizeof(T)*N, cudaMemAttachGlobal );
#else
  ptr = axom::utilities::alloc< T >( N );
#endif
  return ptr;
}

//------------------------------------------------------------------------------
template < typename T >
void deallocate( T*& ptr )
{
  if ( ptr != nullptr )
  {
#if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA)
    cudaFree( ptr );
#else
    axom::utilities::free( ptr );
#endif
  }
}

//------------------------------------------------------------------------------
template< typename execution_policy >
void raja_basic_usage_test( )
{
  constexpr int N = 100;
  int* a = allocate< int >( N );
  int* b = allocate< int >( N );
  int* c = allocate< int >( N );

  // initialize
  RAJA::forall< execution_policy >( RAJA::RangeSegment(0,N), [=](int i) {
      a[ i ] = b[ i ] = 1;
      c[ i ] = 0;
    } );

  // add vectors
  RAJA::forall< execution_policy >( RAJA::RangeSegment(0,N), [=](int i) {
      c[ i ] = a[ i ] + b[ i ];
    } );

  // check result in serial
  for ( int i=0 ; i < N ; ++i )
  {
    EXPECT_EQ( c[ i ], 2 );
  }

  deallocate( a );
  deallocate( b );
  deallocate( c );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST( raja_smoke, basic_use )
{
  std::cout << "Testing RAJA Sequential execution" << std::endl;
  raja_basic_usage_test< RAJA::seq_exec >( );

#if defined(AXOM_USE_OPENMP) && defined(RAJA_ENABLE_OPENMP)
  std::cout << "Testing RAJA OpenMP(CPU) execution" << std::endl;
  raja_basic_usage_test< RAJA::omp_parallel_for_exec >( );
#endif

#if defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA)
  std::cout << "Testing RAJA CUDA execution" << std::endl;
  constexpr int BLOCKSIZE = 256;
  raja_basic_usage_test< RAJA::cuda_exec< BLOCKSIZE > >( );
#endif
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS( );
  return( result );
}
