// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


// Axom primitives
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"

// Axom operations
#include "axom/primal/operators/intersect.hpp"

#include "axom/slic.hpp"

// Testing memory_management
#include "axom/core/memory_management.hpp"

// C++ headers
#include <cmath> // do we need this?
#include <iostream>
#include <fstream>
#include <vector>

#include "fmt/fmt.hpp"

// RAJA 
#include "RAJA/RAJA.hpp"

/*
  CUDA_BLOCK_SIZE - specifies the number of threads in a CUDA thread block
*/
#if defined(RAJA_ENABLE_CUDA)
const int CUDA_BLOCK_SIZE = 256;
#endif

const int in3D = 3;


// primitives represented by doubles in 3D
using PointType = axom::primal::Point<double, in3D>;
using TriangleType = axom::primal::Triangle<double, in3D>;

void makeNonIntersectingTriangles(std::vector<TriangleType> & tris)
{
  PointType p[5];
  p[0] = PointType::make_point(0.0,  1.0, 0.0);
  p[1] = PointType::make_point(-1.0,  0.0, 0.0);
  p[2] = PointType::make_point(1.0,  0.0, 0.0);
  p[3] = PointType::make_point(0.0, 0.0, -1.0);
  p[4] = PointType::make_point(0.0,  0.0, 1.0);

  TriangleType t0(p[0], p[1], p[3]);    tris.push_back(t0);
  TriangleType t1(p[0], p[1], p[4]);    tris.push_back(t1);
  TriangleType t2(p[0], p[2], p[3]);    tris.push_back(t2);
  TriangleType t3(p[0], p[2], p[4]);    tris.push_back(t3);
}

void makeIntersectingTriangles(std::vector<TriangleType> & tris)
{
  PointType p[5];
  p[0] = PointType::make_point(0.0,  1.0, 0.0);
  p[1] = PointType::make_point(-1.0,  0.0, 0.0);
  p[2] = PointType::make_point(1.0,  0.0, 0.0);
  p[3] = PointType::make_point(0.0, 0.0, -1.0);
  p[4] = PointType::make_point(0.0,  0.0, 1.0);

  TriangleType t0(p[0], p[1], p[2]);    tris.push_back(t0);
  TriangleType t1(p[0], p[3], p[4]);    tris.push_back(t1);
  TriangleType t2(p[1], p[2], p[3]);    tris.push_back(t2);
  TriangleType t3(p[1], p[2], p[4]);    tris.push_back(t3);
  TriangleType t4(p[1], p[3], p[4]);    tris.push_back(t4);
  TriangleType t5(p[2], p[3], p[4]);    tris.push_back(t5);
}


void printPairs(std::string title,
                std::vector< std::pair<int, int> > & clashes)
{
  int ccount = clashes.size();
  std::cout << ccount << title << std::endl;
  for (int i = 0 ; i < ccount ; ++i)
  {
    std::cout << clashes[i].first << "   " << clashes[i].second << std::endl;
  }
}

// _naive_triintersect_start
void findTriIntersectionsNaively(
  std::vector<TriangleType> & tris,
  std::vector< std::pair<int, int> > & clashes
  )
{
  int tcount = tris.size();

  for (int i = 0 ; i < tcount ; ++i)
  {
    TriangleType & t1 = tris[i];
    for (int j = i + 1 ; j < tcount ; ++j)
    {
      TriangleType & t2 = tris[j];
      if (intersect(t1, t2))
      {
        clashes.push_back(std::make_pair(i, j));
      }
    }
  }
}
// _naive_triintersect_end

void driveUniformGrid()
{
  // 0 intersections
  std::vector<TriangleType> trisNonIntersect;
  makeNonIntersectingTriangles(trisNonIntersect);

  // 5 intersections
  std::vector<TriangleType> trisIntersect;
  makeIntersectingTriangles(trisIntersect);

  // Stores the clashes (result; needs to be allocated another way)
  std::vector< std::pair<int, int> > naiveclashesNonIntersect;
  findTriIntersectionsNaively(trisNonIntersect, naiveclashesNonIntersect);

  std::vector< std::pair<int, int> > naiveclashesIntersect;
  findTriIntersectionsNaively(trisIntersect, naiveclashesIntersect);

  std::cout << "----- driveUniformGrid -----" << std::endl;
  printPairs(" clashes found by naive algorithm for non-intersecting triangles (0 is expected):", naiveclashesNonIntersect);
  printPairs(" clashes found by naive algorithm for intersecting triangles (5 is expected):", naiveclashesIntersect);
  std::cout << std::endl;
}


void try_axom_umpire(umpire::Allocator allocator)
{
  //Testing memory_management with primitive data type
  int * a = axom::allocate <int> (10, allocator);
  int * b = axom::allocate <int> (10, allocator);
  int * c = axom::allocate <int> (10, allocator);

  for (int i = 0; i < 10; ++i) {
    b[i] = 10;
    a[i] = 1;
    c[i] = 0;
  }

  using EXEC_POL = RAJA::seq_exec;
  //using EXEC_POL = RAJA::cuda_exec<CUDA_BLOCK_SIZE>;

  //add
  RAJA::forall< EXEC_POL >(RAJA::RangeSegment(0, 10), [=] AXOM_HOST_DEVICE (int i) { 
  //RAJA::forall< EXEC_POL >(RAJA::RangeSegment(0, 10), [=] (int i) { 
    c[i] = a[i] + b[i]; 
  });    

  //Print
  std::cout << "----- try_axom_umpire -----" << std::endl;
  for (int i = 0; i < 10; ++i) {
    std::cout << c[i] << std::endl;
  }
  std::cout << std::endl;

  axom::deallocate (a);
  axom::deallocate (b);
  axom::deallocate (c);

}


void try_RAJA_assignment(umpire::Allocator allocator)
{
  // Testing data type
  TriangleType * trie = axom::allocate <TriangleType> (5, allocator);
  PointType * p = axom::allocate <PointType> (5, allocator);
  
  p[0] = PointType::make_point(0.0,  1.0, 0.0);
  p[1] = PointType::make_point(-1.0,  0.0, 0.0);
  p[2] = PointType::make_point(1.0,  0.0, 0.0);
  p[3] = PointType::make_point(0.0, 0.0, -1.0);
  p[4] = PointType::make_point(0.0,  0.0, 1.0);

  TriangleType t0(p[0], p[1], p[2]); 
  TriangleType t1(p[0], p[3], p[4]); 
  TriangleType t2(p[1], p[2], p[3]);
  TriangleType t3(p[1], p[2], p[4]); 
  TriangleType t4(p[1], p[3], p[4]); 

  trie[0] = t0;
  trie[1] = t1;
  trie[2] = t2;
  trie[3] = t3;
  trie[4] = t4;

  using EXEC_POL = RAJA::seq_exec;

  // Attempting to add modification...
  // Works for seq_exec
  RAJA::forall< EXEC_POL >(RAJA::RangeSegment(0, 5), [=] (int i) { 
    //PointType myp = PointType::make_point(10.0,  10.0, 10.0);
    //trie[i][0] = myp;
    trie[i][0] = p[3];
  }); 

  std::cout << "----- try_RAJA_assignment -----" << std::endl;
  std::cout << "Trying assignment in RAJA loop (1st point should all be the same)" << std::endl;
  std::cout << trie[0] << std::endl;
  std::cout << trie[1] << std::endl;
  std::cout << trie[2] << std::endl;
  std::cout << trie[3] << std::endl;
  std::cout << trie[4] << std::endl;
  std::cout << std::endl;

  axom::deallocate(trie);
  axom::deallocate(p);
}



void try_intersect(umpire::Allocator allocator)
{
  //Setup data structures
  TriangleType * tri = axom::allocate <TriangleType> (6, allocator);
  PointType * p = axom::allocate <PointType> (5, allocator);
  
  p[0] = PointType::make_point(0.0,  1.0, 0.0);
  p[1] = PointType::make_point(-1.0,  0.0, 0.0);
  p[2] = PointType::make_point(1.0,  0.0, 0.0);
  p[3] = PointType::make_point(0.0, 0.0, -1.0);
  p[4] = PointType::make_point(0.0,  0.0, 1.0);

  TriangleType t0(p[0], p[1], p[2]); 
  TriangleType t1(p[0], p[3], p[4]); 
  TriangleType t2(p[1], p[2], p[3]);
  TriangleType t3(p[1], p[2], p[4]); 
  TriangleType t4(p[1], p[3], p[4]); 
  TriangleType t5(p[2], p[3], p[4]);

  tri[0] = t0;
  tri[1] = t1;
  tri[2] = t2;
  tri[3] = t3;
  tri[4] = t4;
  tri[5] = t5;

  RAJA::RangeSegment row_range(0, 6);
  RAJA::RangeSegment col_range(0, 6);

  // Sequential Policies
  // using REDUCE_POL = RAJA::seq_reduce;
  // using ATOMIC_POL= RAJA::seq_atomic;
  // using KERNEL_EXEC_POL = 
  //   RAJA::KernelPolicy<
  //     RAJA::statement::For<1, RAJA::seq_exec,
  //       RAJA::statement::For<0, RAJA::seq_exec,
  //         RAJA::statement::Lambda<0>
  //       >
  //     >
  //   >;
  

  // OpenMP Policies
  using REDUCE_POL = RAJA::omp_reduce;
  using ATOMIC_POL= RAJA::omp_atomic;
  using KERNEL_EXEC_POL = 
    RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::omp_parallel_for_exec,
        RAJA::statement::For<0, RAJA::loop_exec,
          RAJA::statement::Lambda<0>
        >
      >
    >;
  

  // CUDA Policies
  // using REDUCE_POL = RAJA::cuda_reduce;
  // using ATOMIC_POL= RAJA::cuda_atomic;
  // using KERNEL_EXEC_POL = 
  //   RAJA::KernelPolicy<
  //     RAJA::statement::CudaKernelFixed<CUDA_BLOCK_SIZE,
  //       RAJA::statement::For<1, RAJA::cuda_block_x_loop,
  //         RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
  //           RAJA::statement::Lambda<0>
  //         >
  //       >
  //     >
  //   >;


  RAJA::ReduceSum< REDUCE_POL, int > numIntersect(0);  

  // RAJA loop to find size to initialize (number of intersections)
    // NOte : reduction more optimal
  RAJA::kernel<KERNEL_EXEC_POL>( RAJA::make_tuple(col_range, row_range),
    [=] AXOM_HOST_DEVICE (int col, int row) {
    if (row > col)
    {
      if (intersect (tri[row], tri[col]))
      {
        numIntersect += 1;
      }
    }
  });
  
  // Print size
  std::cout << "----- try_intersect -----" << std::endl;
  std::cout << "Number of intersections is " << numIntersect.get() << " (Expected was 5)" << std::endl;

  // Allocation to hold intersection pairs and counter to know where to store
  int * intersections = axom::allocate <int> (numIntersect.get() * 2, allocator);
  int * counter = axom::allocate <int> (1, allocator);
  
  counter[0] = 0;

  // RAJA loop to populate with intersections
  RAJA::kernel<KERNEL_EXEC_POL>( RAJA::make_tuple(col_range, row_range),
    [=] AXOM_HOST_DEVICE (int col, int row) {
    if (row > col)
    {
      if (intersect (tri[row], tri[col]))
      {
        auto idx = RAJA::atomicAdd<ATOMIC_POL>(counter, 2);
        intersections[idx] = row;
        intersections[idx+1] = col;

      }
    }
  });

  // Print intersection pairs
  std::cout << "Intersection pairs were: " << std::endl;
  for (int i = 0 ; i < numIntersect.get() * 2 ; i = i + 2)
  {
    std::cout << intersections[i] << "   " << intersections[i + 1] << std::endl;
  }

  // Deallocate 
  axom::deallocate(tri);
  axom::deallocate(p);
  axom::deallocate(intersections);
  axom::deallocate(counter);
}



int main(int argc, char** argv)
{

   #ifdef AXOM_USE_UMPIRE
    std::cout << "UMPIRE IS USED!!!" << std::endl;
  #endif

  #if defined(RAJA_ENABLE_OPENMP)
    std::cout << "OPENMP IS ENABLED!!!" << std::endl;

    umpire::Allocator allocator= axom::getAllocator( axom::getResourceAllocatorID (umpire::resource::Host));
  #endif

  #if defined(RAJA_ENABLE_CUDA)
    std::cout << "CUDA IS ENABLED!!!" << std::endl;
    umpire::Allocator allocator= axom::getAllocator( axom::getResourceAllocatorID (umpire::resource::Unified));
  #endif

  // Deal with unused variables
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);

  

  driveUniformGrid();
  try_axom_umpire(allocator);
  try_RAJA_assignment(allocator);
  try_intersect(allocator);


  return 0;
}
