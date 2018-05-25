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

/*!
 * \file
 *
 * \brief Consists of utility functions to facilitate in test development.
 */
// Axom includes
#include "axom/core/Macros.hpp"

// Mint includes
#include "axom/mint/config.hpp"

#include "axom/mint/mesh/CurvilinearMesh.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/MeshTypes.hpp"
#include "axom/mint/mesh/ParticleMesh.hpp"
#include "axom/mint/mesh/RectilinearMesh.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"

// Slic includes
#include "axom/slic.hpp"

// namespace aliases
namespace mint = axom::mint;

//------------------------------------------------------------------------------
template < int MeshType >
struct mesh_type_name
{
  static constexpr char* name() { return (char*)"[UNDEFINED]"; };
};

//------------------------------------------------------------------------------
template < >
struct mesh_type_name< mint::STRUCTURED_UNIFORM_MESH >
{
  static constexpr char* name() { return (char*)"STRUCTURED_UNIFORM_MESH"; };
};

//------------------------------------------------------------------------------
template < >
struct mesh_type_name< mint::STRUCTURED_CURVILINEAR_MESH >
{
  static constexpr char* name() { return (char*)"STRUCTURED_CURVILINEAR_MESH";};
};

//------------------------------------------------------------------------------
template < >
struct mesh_type_name< mint::STRUCTURED_RECTILINEAR_MESH >
{
  static constexpr char* name() { return (char*)"STRUCTURED_RECTILINEAR_MESH";};
};

//------------------------------------------------------------------------------
template < >
struct mesh_type_name< mint::UNSTRUCTURED_MESH >
{
  static constexpr char* name() { return (char*)"UNSTRUCTURED_MESH";};
};

//------------------------------------------------------------------------------
template < >
struct mesh_type_name< mint::PARTICLE_MESH>
{
  static constexpr char* name() { return (char*)"PARTICLE_MESH";};
};

/*!
 * \brief Creates a mesh of type MeshType equivalent to the given uniform mesh.
 *
 * \param [in]  uniform_mesh the specified uniform mesh.
 * \param [out] output mesh pointer to the output mesh, created by this method.
 *
 * \tparam MeshType the output mesh type
 *
 * \note This method is specialized for each mesh type.
 *
 * \note The caller must deallocate the returned mesh object.
 */
/// @{
template < int MeshType=mint::STRUCTURED_UNIFORM_MESH >
void create_mesh( const mint::UniformMesh* uniform_mesh,
                   mint::Mesh*& output_mesh )
{
  SLIC_ASSERT( uniform_mesh != nullptr );
  SLIC_ASSERT( output_mesh == nullptr );

  const int dimension = uniform_mesh->getDimension();
  mint::int64 ext[ 6 ];

  for ( int i=0; i < dimension; ++i )
  {
    ext[ i*2   ] = 0;
    ext[ i*2+1 ] = uniform_mesh->getNumberOfNodesAlongDim( i )-1;
  }

  output_mesh = new mint::UniformMesh( uniform_mesh->getDimension(),
                                       uniform_mesh->getOrigin(),
                                       uniform_mesh->getSpacing(),
                                       ext );
}

//------------------------------------------------------------------------------
template < >
void create_mesh< mint::STRUCTURED_CURVILINEAR_MESH >(
    const mint::UniformMesh* uniform_mesh, mint::Mesh*& output_mesh )
{
  SLIC_ASSERT( uniform_mesh != nullptr );
  SLIC_ASSERT( output_mesh == nullptr );

  const int dimension         = uniform_mesh->getDimension();
  mint::IndexType node_dims[] = { -1, -1, -1 };

  for ( int i=0; i < dimension; ++i )
  {
    node_dims[ i ] = uniform_mesh->getNumberOfNodesAlongDim( i );
  }

  output_mesh = new mint::CurvilinearMesh( node_dims[ mint::I_DIRECTION ],
                                           node_dims[ mint::J_DIRECTION ],
                                           node_dims[ mint::K_DIRECTION ] );

  const mint::IndexType numNodes = uniform_mesh->getNumberOfNodes();
  for ( mint::IndexType inode=0; inode < numNodes; ++inode )
  {
    double pt[ 3 ];
    uniform_mesh->getNode( inode, pt );

    for ( int idim=0; idim < dimension; ++idim )
    {
      double* coord_array = output_mesh->getCoordinateArray( idim );
      SLIC_ASSERT( coord_array != nullptr );

      coord_array[ inode ] = pt[ idim ];
    }

  } // END for all nodes
}

//------------------------------------------------------------------------------
template < >
void create_mesh< mint::STRUCTURED_RECTILINEAR_MESH >(
    const mint::UniformMesh* uniform_mesh, mint::Mesh*& output_mesh )
{
  SLIC_ASSERT( uniform_mesh != nullptr );
  SLIC_ASSERT( output_mesh == nullptr );

  const int dimension         = uniform_mesh->getDimension();
  SLIC_ASSERT( dimension >= 1 );

  mint::IndexType node_dims[] = { -1, -1, -1 };
  for ( int i=0; i < dimension; ++i )
  {
    node_dims[ i ] = uniform_mesh->getNumberOfNodesAlongDim( i );
  }

  output_mesh = new mint::RectilinearMesh( node_dims[ mint::I_DIRECTION ],
                                           node_dims[ mint::J_DIRECTION ],
                                           node_dims[ mint::K_DIRECTION ]  );

  mint::IndexType Ni =
      uniform_mesh->getNumberOfNodesAlongDim( mint::I_DIRECTION );
  double* x = output_mesh->getCoordinateArray( mint::X_COORDINATE );
  SLIC_ASSERT( x != nullptr );
  for ( mint::IndexType i=0; i < Ni; ++i )
  {
    x[ i ] = uniform_mesh->evaluateCoordinate( i, mint::I_DIRECTION );
  } // END for all i

  if ( dimension >= 2 )
  {
    // fill y
    double* y = output_mesh->getCoordinateArray( mint::Y_COORDINATE );
    SLIC_ASSERT( y != nullptr );

    mint::IndexType Nj =
        uniform_mesh->getNumberOfNodesAlongDim( mint::J_DIRECTION );
    for ( mint::IndexType j=0; j < Nj; ++j )
    {
      y[ j ] = uniform_mesh->evaluateCoordinate( j, mint::J_DIRECTION );
    } // END for all j

  }

  if ( dimension == 3 )
  {
    // fill z
    double* z = output_mesh->getCoordinateArray( mint::Z_COORDINATE );
    SLIC_ASSERT( z != nullptr );

    const mint::IndexType Nk =
        uniform_mesh->getNumberOfNodesAlongDim( mint::K_DIRECTION );
    for ( mint::IndexType k=0; k < Nk; ++k )
    {
      z[ k ] = uniform_mesh->evaluateCoordinate( k, mint::K_DIRECTION );
    } // END for all k

  } // END if 3-D

}

//------------------------------------------------------------------------------
template < >
void create_mesh< mint::PARTICLE_MESH >(
    const mint::UniformMesh* uniform_mesh, mint::Mesh*& output_mesh )
{
  SLIC_ASSERT( uniform_mesh != nullptr );
  SLIC_ASSERT( output_mesh == nullptr );

  const int dimension            = uniform_mesh->getDimension();
  const mint::IndexType numNodes = uniform_mesh->getNumberOfNodes();

  output_mesh = new mint::ParticleMesh( dimension, numNodes );

  for ( mint::IndexType inode=0; inode < numNodes; ++inode )
  {
    double node[ 3 ];
    uniform_mesh->getNode( inode, node );

    for ( int idim=0; idim < dimension; ++idim )
    {
      double* coord = output_mesh->getCoordinateArray( idim );
      SLIC_ASSERT( coord != nullptr );

      coord[ inode ] = node[ idim ];
    }

  } // END for all nodes

}

//------------------------------------------------------------------------------
template < >
void create_mesh< mint::UNSTRUCTURED_MESH >(
    const mint::UniformMesh* uniform_mesh, mint::Mesh*& output_mesh )
{
  SLIC_ASSERT( uniform_mesh != nullptr );
  SLIC_ASSERT( output_mesh == nullptr );

  const int dimension            = uniform_mesh->getDimension();
  const mint::IndexType numNodes = uniform_mesh->getNumberOfNodes();
  const mint::IndexType numCells = uniform_mesh->getNumberOfCells();

  using UnstructuredMeshType = mint::UnstructuredMesh< mint::SINGLE_SHAPE >;
  mint::CellType cell_type   = ( dimension==3 )? mint::HEX :
                                 ( (dimension==2)? mint::QUAD : mint::SEGMENT );

  output_mesh = new UnstructuredMeshType( dimension, cell_type,
                                          numNodes, numCells );

  UnstructuredMeshType* m = static_cast< UnstructuredMeshType* >( output_mesh );

  // append nodes
  for ( mint::IndexType inode=0; inode < numNodes; ++inode )
  {
    double node[ 3 ];
    uniform_mesh->getNode( inode, node );
    m->appendNodes( node );
  }

  // append cells
  for ( mint::IndexType icell=0; icell < numCells; ++icell )
  {
    mint::IndexType cell[ 8 ];
    uniform_mesh->getCell( icell, cell );
    m->appendCell( cell, cell_type );
  }

}
/// @}
