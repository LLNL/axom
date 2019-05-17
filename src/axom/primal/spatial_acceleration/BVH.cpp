// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/Types.hpp"

#include "axom/primal/spatial_acceleration/BVH.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/aabb.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/vec.hpp"
#include "axom/primal/spatial_acceleration/linear_bvh/policies.hpp"


// RAJA includes
#include "RAJA/RAJA.hpp"

#include "axom/slic/interface/slic.hpp" // for SLIC macros

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

// C/C++ includes
#include <fstream> // for std::ofstream
#include <sstream> // for std::ostringstream

namespace axom
{
namespace primal
{

namespace
{
void find_candidates( IndexType* offsets,
                       IndexType*& candidates,
                       IndexType numPts,
                       const double* x,
                       const double* y,
                       const double* z,
                       bvh::BVH<float32,3> &bvh)
{


  IndexType *candidate_counts = axom::allocate<IndexType>(numPts);
  const bvh::Vec<float32, 4>  *inner_nodes = bvh.m_inner_nodes;
  const int32 *leaf_nodes = bvh.m_leaf_nodes;
  RAJA::forall<bvh::raja_for_policy>(RAJA::RangeSegment(0, numPts), AXOM_LAMBDA (IndexType i)
  {
    int32 count = 0;

    bvh::Vec<double, 3> point;
    point[0] = x[i];
    point[1] = y[i];
    point[2] = z[i];

    int32 current_node = 0;
    int32 todo[64];
    int32 stackptr = 0;

    constexpr int32 barrier = -2000000000;
    todo[stackptr] = barrier;
    while (current_node != barrier)
    {
      if (current_node > -1)
      {
        // inner node
        // TODO: implement const get vec4
        //const bvh::Vec<float32, 4> first4  = const_get_vec4f(&inner_nodes[current_node + 0]);
        //const bvh::Vec<float32, 4> second4 = const_get_vec4f(&inner_nodes[current_node + 1]);
        //const bvh::Vec<float32, 4> third4  = const_get_vec4f(&inner_nodes[current_node + 2]);
        const bvh::Vec<float32, 4> first4  = inner_nodes[current_node + 0];
        const bvh::Vec<float32, 4> second4 = inner_nodes[current_node + 1];
        const bvh::Vec<float32, 4> third4  = inner_nodes[current_node + 2];

        bool in_left = true;
        if(point[0]  < first4[0]) in_left = false; // left x min
        if(point[1]  < first4[1]) in_left = false; // y min
        if(point[2]  < first4[2]) in_left = false; // z min

        if(point[0]  > first4[3])  in_left = false; // x max
        if(point[1]  > second4[0]) in_left = false;
        if(point[2]  > second4[1]) in_left = false;

        bool in_right = true;
        if(point[0]  < second4[2]) in_right = false;
        if(point[1]  < second4[3]) in_right = false;
        if(point[2]  < third4[0])  in_right = false;

        if(point[0]  > third4[1]) in_right = false;
        if(point[1]  > third4[2]) in_right = false;
        if(point[2]  > third4[3]) in_right = false;

        if (!in_left && !in_right)
        {
          // pop the stack and continue
          current_node = todo[stackptr];
          stackptr--;
        }
        else
        {
          // TODO: implement const get vec4
          //bvh::Vec<float32, 4> children = const_get_vec4f(&inner_ptr[current_node + 3]);
          bvh::Vec<float32, 4> children = inner_nodes[current_node + 3];
          int32 l_child;
          constexpr int32 isize = sizeof(int32);
          // memcpy the int bits hidden in the floats
          memcpy(&l_child, &children[0], isize);
          int32 r_child;
          memcpy(&r_child, &children[1], isize);

          current_node = (in_left) ? l_child : r_child;

          if (in_left && in_right)
          {
            stackptr++;
            todo[stackptr] = r_child;
            // TODO: if we are in both children we could
            // go down the "closer" first by perhaps the distance
            // from the point to the center of the aabb
          }
        }
      }
      else
      {
        // leaf node
        count++;
        current_node = todo[stackptr];
        stackptr--;
      }
    } // while
    candidate_counts[i] = count;
//#endif
  });

  RAJA::exclusive_scan<bvh::raja_for_policy>(candidate_counts,
                                             candidate_counts + numPts,
                                             offsets,
                                             RAJA::operators::plus<IndexType>{});

  // TODO: this will segault with raw(unmanaged) cuda pointers
  IndexType total_candidates = offsets[numPts - 1] + candidate_counts[numPts - 1];

  candidates = axom::allocate< IndexType >( total_candidates);

  //candidates =
  RAJA::forall<bvh::raja_for_policy>(RAJA::RangeSegment(0, numPts), AXOM_LAMBDA (IndexType i)
  {
    int32 offset = offsets[i];

    bvh::Vec<double, 3> point;
    point[0] = x[i];
    point[1] = y[i];
    point[2] = z[i];

    int32 current_node = 0;
    int32 todo[64];
    int32 stackptr = 0;

    constexpr int32 barrier = -2000000000;
    todo[stackptr] = barrier;
    while (current_node != barrier)
    {
      if (current_node > -1)
      {
        // inner node
        // TODO: implement const get vec4
        //const bvh::Vec<float32, 4> first4  = const_get_vec4f(&inner_nodes[current_node + 0]);
        //const bvh::Vec<float32, 4> second4 = const_get_vec4f(&inner_nodes[current_node + 1]);
        //const bvh::Vec<float32, 4> third4  = const_get_vec4f(&inner_nodes[current_node + 2]);
        const bvh::Vec<float32, 4> first4  = inner_nodes[current_node + 0];
        const bvh::Vec<float32, 4> second4 = inner_nodes[current_node + 1];
        const bvh::Vec<float32, 4> third4  = inner_nodes[current_node + 2];

        bool in_left = true;
        if(point[0]  < first4[0]) in_left = false;
        if(point[1]  < first4[1]) in_left = false;
        if(point[2]  < first4[2]) in_left = false;

        if(point[0]  > first4[3])  in_left = false;
        if(point[1]  > second4[0]) in_left = false;
        if(point[2]  > second4[1]) in_left = false;

        bool in_right = true;
        if(point[0]  < second4[2]) in_right = false;
        if(point[1]  < second4[3]) in_right = false;
        if(point[2]  < third4[0])  in_right = false;

        if(point[0]  > third4[1]) in_right = false;
        if(point[1]  > third4[2]) in_right = false;
        if(point[2]  > third4[3]) in_right = false;

        if (!in_left && !in_right)
        {
          // pop the stack and continue
          current_node = todo[stackptr];
          stackptr--;
        }
        else
        {
          // TODO: implement const get vec4
          //bvh::Vec<float32, 4> children = const_get_vec4f(&inner_ptr[current_node + 3]);
          bvh::Vec<float32, 4> children = inner_nodes[current_node + 3];
          int32 l_child;
          constexpr int32 isize = sizeof(int32);
          // memcpy the int bits hidden in the floats
          memcpy(&l_child, &children[0], isize);
          int32 r_child;
          memcpy(&r_child, &children[1], isize);

          current_node = (in_left) ? l_child : r_child;

          if (in_left && in_right)
          {
            stackptr++;
            todo[stackptr] = r_child;
            // TODO: if we are in both children we could
            // go down the "closer" first by perhaps the distance
            // from the point to the center of the aabb
          }
        }
      }
      else
      {
        current_node = -current_node - 1; //swap the neg address
        candidates[offset] = leaf_nodes[current_node];
        offset++;
        current_node = todo[stackptr];
        stackptr--;
      }
    } // while
//#endif
  });

  axom::deallocate(candidate_counts);
}

//------------------------------------------------------------------------------
void write_box( const float32* min,
               const float32* max,
               int& numPoints,
               int& numBins,
               std::ostringstream& nodes,
               std::ostringstream& cells )
{
  nodes << min[0] << " " << min[1] << " " << min[2] << std::endl;
  nodes << max[0] << " " << min[1] << " " << min[2] << std::endl;
  nodes << max[0] << " " << max[1] << " " << min[2] << std::endl;
  nodes << min[0] << " " << max[1] << " " << min[2] << std::endl;

  nodes << min[0] << " " << min[1] << " " << max[2] << std::endl;
  nodes << max[0] << " " << min[1] << " " << max[2] << std::endl;
  nodes << max[0] << " " << max[1] << " " << max[2] << std::endl;
  nodes << min[0] << " " << max[1] << " " << max[2] << std::endl;

  int offset = numPoints;
  cells << "8 ";
  for ( int i=0; i < 8; ++i )
  {
    cells << offset + i << " ";
  }
  cells << "\n";
  ++numBins;

  numPoints += 8;
}

//------------------------------------------------------------------------------
void write_bvh_bins_recursive( bvh::Vec<float32, 4>* inner_nodes,
                               int32 current_node, int level,
                               int& numPoints,
                               int& numBins,
                               std::ostringstream& nodes,
                               std::ostringstream& cells,
                               std::ostringstream& levels )
{
  // STEP 0: get the Left and Right boxes in the flat layout
  const bvh::Vec<float32, 4> first4  = inner_nodes[current_node + 0];
  const bvh::Vec<float32, 4> second4 = inner_nodes[current_node + 1];
  const bvh::Vec<float32, 4> third4  = inner_nodes[current_node + 2];

  // STEP 1: get the right and left child
  int32 l_child;
  int32 r_child;
  constexpr int32 isize = sizeof(int32);
  bvh::Vec<float32, 4> children = inner_nodes[current_node + 3];
  // memcpy the int bits hidden in the floats
  memcpy(&l_child, &children[0], isize);
  memcpy(&r_child, &children[1], isize);

  // STEP 2: check left
  if ( l_child > -1 )
  {
    float32 min[ 3 ] = { first4[0], first4[1], first4[2]   };
    float32 max[ 3 ] = { first4[3], second4[0], second4[1] };

    write_box( min, max, numPoints, numBins, nodes, cells );
    levels << level << std::endl;

    write_bvh_bins_recursive( inner_nodes, l_child, level+1,
                              numPoints, numBins,
                              nodes,
                              cells,
                              levels );
  }

  // STEP 3: check right
  if ( r_child > -1 )
  {
    float32 min[ 3 ] = { second4[2], second4[3], third4[0] };
    float32 max[ 3 ] = { third4[1], third4[2], third4[3] };

    write_box( min, max, numPoints, numBins, nodes, cells );
    levels << level << std::endl;

    write_bvh_bins_recursive( inner_nodes, r_child, level+1,
                              numPoints, numBins,
                              nodes,
                              cells,
                              levels );
  }

}

//------------------------------------------------------------------------------
void write_bvh_root( const bvh::AABB<float32,3>& root,
                     int& numPoints,
                     int& numBins,
                     std::ostringstream& nodes,
                     std::ostringstream& cells,
                     std::ostringstream& levels  )
{
  float32 min[ 3 ];
  float32 max[ 3 ];

  min[ 0 ] = root.m_x.min();
  min[ 1 ] = root.m_y.min();
  min[ 2 ] = root.m_z.min();

  max[ 0 ] = root.m_x.max();
  max[ 1 ] = root.m_y.max();
  max[ 2 ] = root.m_z.max();

  int32 offset  = numPoints;
  write_box( min, max, numPoints, numBins, nodes, cells );
  levels << "0\n";
}

} /* namespace anonymous*/

//------------------------------------------------------------------------------
// BVH IMPLEMENTATION
//------------------------------------------------------------------------------

BVH::BVH( int dimension, const double* boxes, IndexType numItems ) :
    m_dimension( dimension ),
    m_numItems( numItems ),
    m_boxes( boxes )
{
  SLIC_ASSERT( m_dimension >= 1 && m_dimension <= 3 );
  SLIC_ASSERT( m_boxes != nullptr );
  SLIC_ASSERT( m_numItems >= 1 );

  // TODO: implement this
}

//------------------------------------------------------------------------------
BVH::~BVH()
{
  m_bvh.free();
}

//------------------------------------------------------------------------------
int BVH::build( )
{

  bvh::LinearBVHBuilder builder;
  m_bvh = builder.construct<float32,3>( static_cast<float32*>(m_boxes), m_numItems);
  std::cout<<"BOUNDS "<<m_bvh.m_bounds<<"\n";
  return BVH_BUILD_OK;
}

//------------------------------------------------------------------------------
void BVH::find( IndexType*& candidates,
                IndexType& numCandidates,
                double x,
                double y,
                double z )
{
  SLIC_ASSERT( candidates == nullptr );

  numCandidates = 0;
  //find_candidates(offsets, &candidates, 1, &x, &y, &z, m_bvh);
  // TODO: implement this
}

//------------------------------------------------------------------------------
void BVH::find( IndexType* offsets,
                IndexType*& candidates,
                IndexType numPts,
                const double* x,
                const double* y,
                const double* z )
{
  find_candidates(offsets, candidates, numPts, x, y, z, m_bvh);
  // TODO: implement this
}

//------------------------------------------------------------------------------
void BVH::writeVtkFile( const std::string& fileName ) const
{
//  int32 numPoints = 0;
//  int32 numBins   = 0;
//
//  std::ostringstream coordinates;
//  std::ostringstream cells;
//  std::ostringstream levels;
//
//  // STEP 0: Write VTK header
//  std::ofstream ofs;
//  ofs.open( fileName.c_str() );
//  ofs << "# vtk DataFile Version 3.0\n";
//  ofs << " BVHTree \n";
//  ofs << "ASCII\n";
//  ofs << "DATASET UNSTRUCTURED_GRID\n";
//
//  // STEP 1: extract root
//  write_bvh_root( m_bvh.m_bounds, numPoints, numBins,
//                  coordinates, cells, levels );
//  SLIC_ASSERT( numPoints == 8 );
//  SLIC_ASSERT( numBins == 1 );
//
//  // STEP 2: recursively traverse the tree write the data
//  constexpr int32 ROOT = 0;
//  write_bvh_bins_recursive( m_bvh.m_inner_nodes, ROOT, 1, numPoints, numBins,
//                            coordinates,
//                            cells,
//                            levels );
//
//  // STEP 3: write nodal coordinates
//  ofs << "POINTS " << numPoints << " double\n";
//  ofs << coordinates.str() << std::endl;
//
//  // STEP 4: write cell connectivity
//  const int nnodes = (m_dimension==2) ? 4 : 8;
//  ofs << "CELLS " << numBins << " " << numBins*(nnodes+1) << std::endl;
//  ofs << cells.str() << std::endl;
//
//  // STEP 5: write cell types
//  ofs << "CELL_TYPES " << numBins << std::endl;
//  const int cellType = (m_dimension==2) ? 9 : 12;
//  for ( int i=0; i < numBins; ++i )
//  {
//    ofs << cellType << std::endl;
//  } // END for all bins
//
//  // STEP 6: write level information
//  ofs << "CELL_DATA " << numBins << std::endl;
//  ofs << "SCALARS level int\n";
//  ofs << "LOOKUP_TABLE default\n";
//  ofs << levels.str() << std::endl;
//  ofs << std::endl;
//
//  ofs.close();
}

} /* namespace primal */
} /* namespace axom */
