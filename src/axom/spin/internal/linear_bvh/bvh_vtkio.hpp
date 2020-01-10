// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_LINEAR_BVH_VTKIO_HPP_
#define AXOM_SPIN_LINEAR_BVH_VTKIO_HPP_

#include "axom/spin/internal/linear_bvh/vec.hpp"
#include "axom/spin/internal/linear_bvh/aabb.hpp"

namespace axom
{
namespace spin
{
namespace internal
{
namespace linear_bvh
{

//------------------------------------------------------------------------------
template < typename FloatType >
void write_box2d( const FloatType& xmin,
                  const FloatType& ymin,
                  const FloatType& xmax,
                  const FloatType& ymax,
                  int32& numPoints,
                  int32& numBins,
                  std::ostringstream& nodes,
                  std::ostringstream& cells  )
{
  nodes << xmin << " " << ymin << " 0.0\n";
  nodes << xmax << " " << ymin << " 0.0\n";
  nodes << xmax << " " << ymax << " 0.0\n";
  nodes << xmin << " " << ymax << " 0.0\n";

  constexpr int32 NUM_NODES = 4;
  int32 offset = numPoints;
  cells << NUM_NODES << " ";
  for ( int32 i=0; i < NUM_NODES; ++i )
  {
    cells << (offset + i) << " ";
  }
  cells << "\n";

  numBins   += 1;
  numPoints += NUM_NODES;
}

//------------------------------------------------------------------------------
template < typename FloatType >
void write_box3d( const FloatType& xmin,
                  const FloatType& ymin,
                  const FloatType& zmin,
                  const FloatType& xmax,
                  const FloatType& ymax,
                  const FloatType& zmax,
                  int32& numPoints,
                  int32& numBins,
                  std::ostringstream& nodes,
                  std::ostringstream& cells  )
{
  nodes << xmin << " " << ymin << " " << zmin << std::endl;
  nodes << xmax << " " << ymin << " " << zmin << std::endl;
  nodes << xmax << " " << ymax << " " << zmin << std::endl;
  nodes << xmin << " " << ymax << " " << zmin << std::endl;

  nodes << xmin << " " << ymin << " " << zmax << std::endl;
  nodes << xmax << " " << ymin << " " << zmax << std::endl;
  nodes << xmax << " " << ymax << " " << zmax << std::endl;
  nodes << xmin << " " << ymax << " " << zmax << std::endl;

  constexpr int32 NUM_NODES = 8;
  int32 offset = numPoints;
  cells << NUM_NODES << " ";
  for ( int32 i=0; i < NUM_NODES; ++i )
  {
    cells << (offset + i) << " ";
  }
  cells << "\n";

  numBins   += 1;
  numPoints += NUM_NODES;
}

//------------------------------------------------------------------------------
template < typename FloatType, int NDIMS >
void write_leftbox( const Vec< FloatType, 4 >& first,
                    const Vec< FloatType, 4 >& second,
                    int32& numPoints,
                    int32& numBins,
                    std::ostringstream& nodes,
                    std::ostringstream& cells  )
{
  const FloatType& xmin = first[ 0 ];
  const FloatType& ymin = first[ 1 ];

  const FloatType& xmax = first[ 3 ];
  const FloatType& ymax = second[ 0 ];

  if ( NDIMS== 2 )
  {
    write_box2d( xmin, ymin, xmax, ymax,
                 numPoints, numBins, nodes, cells );
  }
  else
  {
    const FloatType& zmin = first[ 2 ];
    const FloatType& zmax = second[ 2 ];

    write_box3d( xmin, ymin, zmin,
                 xmax, ymax, zmax,
                 numPoints, numBins, nodes, cells );
  }

}

//------------------------------------------------------------------------------
template < typename FloatType, int NDIMS >
void write_righbox( const Vec< FloatType, 4 >& second,
                    const Vec< FloatType, 4 >& third,
                    int32& numPoints,
                    int32& numBins,
                    std::ostringstream& nodes,
                    std::ostringstream& cells  )
{
  const FloatType& xmin = second[ 2 ];
  const FloatType& ymin = second[ 3 ];

  const FloatType& xmax = third[ 1 ];
  const FloatType& ymax = third[ 2 ];

  if ( NDIMS== 2 )
  {
    write_box2d( xmin, ymin, xmax, ymax,
                 numPoints, numBins, nodes, cells );
  }
  else
  {
    const FloatType& zmin = third[ 0 ];
    const FloatType& zmax = third[ 3 ];

    write_box3d( xmin, ymin, zmin,
                 xmax, ymax, zmax,
                 numPoints, numBins, nodes, cells );
  }

}

//------------------------------------------------------------------------------
template < typename FloatType >
void write_root( const AABB< FloatType, 2 >& root,
                 int32& numPoints,
                 int32& numBins,
                 std::ostringstream& nodes,
                 std::ostringstream& cells,
                 std::ostringstream& levels )
{

  write_box2d( root.m_x.min(), root.m_y.min(),
               root.m_x.max(), root.m_y.max(),
               numPoints,
               numBins,
               nodes,
               cells  );
  levels << "0\n";

}

//------------------------------------------------------------------------------
template < typename FloatType >
void write_root( const AABB< FloatType, 3 >& root,
                 int32& numPoints,
                 int32& numBins,
                 std::ostringstream& nodes,
                 std::ostringstream& cells,
                 std::ostringstream& levels )
{

  write_box3d( root.m_x.min(), root.m_y.min(), root.m_z.min(),
               root.m_x.max(), root.m_y.max(), root.m_z.max(),
               numPoints,
               numBins,
               nodes,
               cells  );
  levels << "0\n";
}

//------------------------------------------------------------------------------
template < typename FloatType, int NDIMS  >
void write_recursive( Vec< FloatType, 4 >* inner_nodes,
                      int32 current_node,
                      int32 level,
                      int32& numPoints,
                      int32& numBins,
                      std::ostringstream& nodes,
                      std::ostringstream& cells,
                      std::ostringstream& levels )
{
  // STEP 0: get the flat BVH layout
  const Vec< FloatType, 4> first4  = inner_nodes[current_node + 0];
  const Vec< FloatType, 4> second4 = inner_nodes[current_node + 1];
  const Vec< FloatType, 4> third4  = inner_nodes[current_node + 2];

  // STEP 1: extract children information
  int32 l_child;
  int32 r_child;
  constexpr int32 isize = sizeof(int32);
  Vec< FloatType, 4 > children = inner_nodes[current_node + 3];
  memcpy(&l_child,&children[0],isize);
  memcpy(&r_child,&children[1],isize);

  write_leftbox< FloatType, NDIMS >( first4, second4, numPoints, numBins,
                                     nodes, cells );
  levels << level << std::endl;
  write_righbox< FloatType, NDIMS >( second4, third4, numPoints, numBins,
                                     nodes, cells );
  levels << level << std::endl;

  // STEP 2: check left
  if ( l_child > - 1 )
  {
    write_recursive< FloatType, NDIMS >(
        inner_nodes, l_child, level+1, numPoints, numBins,
        nodes, cells, levels );

  }

  // STEP 3: check right
  if ( r_child > -1 )
  {
    write_recursive< FloatType, NDIMS >(
        inner_nodes, r_child, level+1, numPoints, numBins,
        nodes, cells, levels );

  }

}


} /* namespace linear_bvh */
} /* namespace internal */
} /* namespace spin */
} /* namespace axom */

#endif /* AXOM_SPIN_LINEAR_BVH_VTKIO_HPP_ */
