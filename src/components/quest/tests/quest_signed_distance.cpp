/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

#include "quest/BoundingBox.hpp"
#include "quest/CellType.hpp"
#include "quest/HyperSphere.hpp"
#include "quest/MeshType.hpp"
#include "quest/SignedDistance.hpp"
#include "quest/UniformMesh.hpp"
#include "quest/UnstructuredMesh.hpp"

#include "sphere.hpp"

#include "gtest/gtest.h"

typedef meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE > TriangleMesh;
typedef meshtk::UniformMesh UniformMesh;

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace detail
{

/*!
 *******************************************************************************
 * \brief Gets a surface mesh instance for the sphere.
 * \param [in] mesh pointer to the mesh instance.
 * \pre mesh != ATK_NULLPTR
 *******************************************************************************
 */
void getMesh( TriangleMesh* mesh )
{
  SLIC_ASSERT( mesh != ATK_NULLPTR );

  for ( int inode=0; inode < sphere::num_nodes; ++inode ) {
     mesh->insertNode( sphere::nodes[ inode*3   ],
                       sphere::nodes[ inode*3+1 ],
                       sphere::nodes[ inode*3+2 ]  );
  } // END for all nodes

  for ( int icell=0; icell < sphere::num_cells; ++icell ) {
     mesh->insertCell( &sphere::cells[ icell*3 ], meshtk::LINEAR_TRIANGLE,3);
  } // END for all cells

}

/*!
 *******************************************************************************
 * \brief Returns the bounding box of the mesh.
 * \param [in] mesh pointer to the mesh instance.
 * \return bb bounding box of the mesh
 *******************************************************************************
 */
quest::BoundingBox< double,3 > getBounds( const meshtk::Mesh* mesh )
{
  SLIC_ASSERT( mesh != ATK_NULLPTR );

  quest::BoundingBox< double,3 > bb;
  quest::Point< double,3 > pt;

  const int nnodes = mesh->getMeshNumberOfNodes();
  for ( int inode=0; inode < nnodes; ++inode ) {
     mesh->getMeshNode( inode, pt.data() );
     bb.addPoint( pt );
  }

  return( bb );
}

/*!
 *******************************************************************************
 * \brief Generates a uniform mesh surrounding the given triangle mesh.
 * \param [in] mesh pointer to the input mesh.
 * \param [in] umesh pointer to the uniform mesh;
 *******************************************************************************
 */
void getUniformMesh( const TriangleMesh* mesh, UniformMesh*& umesh )
{
  SLIC_ASSERT( mesh != ATK_NULLPTR );
  SLIC_ASSERT( umesh == ATK_NULLPTR );

  quest::BoundingBox< double,3 > bb = getBounds( mesh );
  bb.expand( 2.0 );

  double h[3];
  h[0] = (bb.getMax()[0]-bb.getMin()[0]) / 32;
  h[1] = (bb.getMax()[1]-bb.getMin()[1]) / 32;
  h[2] = (bb.getMax()[2]-bb.getMin()[2]) / 32;

  int ext[6];
  ext[0] = ext[2] = ext[4] = 0;
  ext[1] = ext[3] = ext[5] = 32;

  umesh = new UniformMesh(3,bb.getMin().data(),h,ext);
}

} /* end detail namespace */

//------------------------------------------------------------------------------
TEST( quest_signed_distance, sphere_test )
{
  SLIC_INFO( "Constructing sphere mesh..." );
  TriangleMesh* surface_mesh = new TriangleMesh( 3 );
  detail::getMesh( surface_mesh );

  SLIC_INFO( "Generating uniform mesh..." );
  UniformMesh* umesh = ATK_NULLPTR;
  detail::getUniformMesh( surface_mesh, umesh );

  SLIC_INFO( "Generate BVHTree..." );
  quest::SignedDistance< 3 > signed_distance( surface_mesh, 25, 32 );

  SLIC_INFO( "Compute signed distance..." );
  quest::HyperSphere< double,3 > analytic_sphere( sphere::radius );
  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf   = std::numeric_limits< double >::min( );

  const int nnodes = umesh->getNumberOfNodes();
  for ( int inode=0; inode < nnodes; ++inode ) {

     quest::Point< double,3 > pt;
     umesh->getNode( inode, pt.data() );

     double computed = signed_distance.computeDistance( pt );
     double exact    = analytic_sphere.getSignedDistance( pt.data() );

     // compute error
     double dx    = computed - exact;
     double absdx = std::fabs( dx );

     EXPECT_NEAR( exact, computed, 1.e-2 );

     l1norm += absdx;
     l2norm += dx*dx;

     if ( absdx > linf ) {
       linf = absdx;
     }

  } // END for all nodes

  l2norm = std::sqrt( l2norm );

  SLIC_INFO( "l1 =" << l1norm );
  SLIC_INFO( "l2 =" << l2norm );
  SLIC_INFO( "linf = " << linf );

  EXPECT_NEAR( 20.8774,  l1norm, 1.e-4 );
  EXPECT_NEAR( 0.140352, l2norm, 1.e-4 );
  EXPECT_NEAR( 0.00194587, linf, 1.e-4 );

  delete surface_mesh;
  delete umesh;

  SLIC_INFO( "Done." );
}

//------------------------------------------------------------------------------
int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
