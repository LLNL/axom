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
#include "mint/ParticleMesh.hpp"

#include "axom/Macros.hpp"

#include "mint/config.hpp"          // for compile-time type definitions
#include "mint/MeshCoordinates.hpp" // for mint::MeshCoordinates class
#include "mint/blueprint.hpp"       // for mint::blueprint() functions

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{

inline int dim( const double* AXOM_NOT_USED(x),
                const double* y, const double* z )
{
  return ( ( z != AXOM_NULLPTR ) ? 3 : ( (y != AXOM_NULLPTR ) ? 2 : 1 ) );
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//  PARTICLE MESH IMPLEMENTATION
//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension,
                            IndexType numParticles,
                            IndexType capacity ) :
          ParticleMesh( dimension, 0, 0, numParticles, capacity )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, int blockId, int partId,
                            IndexType numParticles,
                            IndexType capacity ) :
            Mesh( dimension, PARTICLE_MESH, blockId, partId ),
            m_positions( new MeshCoordinates(dimension,numParticles,capacity) )
{
  initialize( );
}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( IndexType numParticles,
                            double* x, double* y, double* z ) :
      Mesh( dim(x,y,z), PARTICLE_MESH, 0, 0 ),
      m_positions( new MeshCoordinates( numParticles, numParticles, x, y, z ) )
{
  initialize( );
}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int blockId, int partId,
                            IndexType numParticles,
                            double* x, double* y, double* z ) :
      Mesh( dim(x,y,z), PARTICLE_MESH, blockId, partId ),
      m_positions( new MeshCoordinates( numParticles, numParticles, x, y, z ) )
{
  initialize( );
}

//------------------------------------------------------------------------------
#ifdef MINT_USE_SIDRE

ParticleMesh::ParticleMesh( sidre::Group* group,
                            const std::string& topo ) :
    Mesh( group, topo ),
    m_positions( new MeshCoordinates( getCoordsetGroup() ) )
{
  SLIC_ERROR_IF( m_type != PARTICLE_MESH,
      "supplied Sidre group does not correspond to a ParticleMesh" );

  initialize( );
}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( sidre::Group* group ) : ParticleMesh( group, "" )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, int blockId, int partId,
                            IndexType numParticles,
                            sidre::Group* group,
                            const std::string& topo,
                            const std::string& coordset,
                            IndexType capacity ) :
    Mesh( dimension, PARTICLE_MESH, blockId, partId, group, topo, coordset )

{
  SLIC_ERROR_IF( !blueprint::validTopologyGroup( getTopologyGroup() ),
                 "invalid topology group!" );

  getTopologyGroup()->getView( "type" )->setString( "particle" );

  m_positions = new MeshCoordinates( getCoordsetGroup(), dimension,
                                     numParticles, capacity );

  initialize( );
}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension,
                            IndexType numParticles,
                            sidre::Group* group,
                            const std::string& topo,
                            const std::string& coordset,
                            IndexType capacity ) :
 ParticleMesh( dimension, 0, 0, numParticles, group, topo, coordset, capacity )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, int blockId, int partId,
                            IndexType numParticles,
                            sidre::Group* group,
                            IndexType capacity ) :
  ParticleMesh(dimension,blockId,partId, numParticles, group, "", "",capacity )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension,
                            IndexType numParticles,
                            sidre::Group* group,
                            IndexType capacity ) :
    ParticleMesh( dimension, 0, 0, numParticles, group, "", "",capacity )
{

}

#endif

//------------------------------------------------------------------------------
void ParticleMesh::initialize( )
{
  SLIC_ASSERT( m_positions != AXOM_NULLPTR );

  m_num_nodes             = m_positions->numNodes();
  m_num_cells             = m_num_nodes;
  m_explicit_coords       = true;
  m_explicit_connectivity = false;
  m_has_mixed_topology    = false;
}

//------------------------------------------------------------------------------
ParticleMesh::~ParticleMesh()
{
  delete m_positions;
  m_positions = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
bool ParticleMesh::checkConsistency()
{
  bool status = ( m_num_cells==m_num_nodes );
  status = status && (m_positions->numNodes()==m_num_nodes);

  const int numFields = m_mesh_fields[ NODE_CENTERED ]->getNumFields();
  for ( int i=0; status && (i < numFields); ++i )
  {
    mint::Field* f = m_mesh_fields[ NODE_CENTERED ]->getField( i );
    status = status && ( f != AXOM_NULLPTR );
    status = status && ( f->getNumTuples()==m_num_nodes );
  }

  return status;
}

} /* namespace mint */
} /* namespace axom */
