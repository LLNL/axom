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

#include "mint/MeshHelpers.hpp"     // for internal helper methods

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, IndexType numParticles,
                            IndexType capacity ) :
  Mesh( dimension, PARTICLE_MESH ),
  m_positions( new MeshCoordinates(dimension, numParticles, capacity) )
{
  initialize( );
}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( IndexType numParticles,
                            double* x, double* y, double* z ) :
  Mesh( internal::dim( x, y, z ), PARTICLE_MESH ),
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
ParticleMesh::ParticleMesh( int dimension, IndexType numParticles,
                            sidre::Group* group,
                            const std::string& topo,
                            const std::string& coordset,
                            IndexType capacity ) :
  Mesh( dimension, PARTICLE_MESH, group, topo, coordset ),
  m_positions( AXOM_NULLPTR )
{
  blueprint::initializeTopologyGroup( m_group, m_topology, m_coordset,
                                      "points" );

  SLIC_ERROR_IF( !blueprint::isValidTopologyGroup( getTopologyGroup() ),
                 "invalid topology group!" );

  m_positions = new MeshCoordinates( getCoordsetGroup(), dimension,
                                     numParticles, capacity );

  initialize( );
}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, IndexType numParticles,
                            sidre::Group* group, IndexType capacity ) :
  ParticleMesh(dimension, numParticles, group, "", "",capacity )
{}

#endif

//------------------------------------------------------------------------------
void ParticleMesh::initialize( )
{
  SLIC_ASSERT( m_positions != AXOM_NULLPTR );

  m_mesh_fields[ NODE_CENTERED ]->reserve( getNodeCapacity() );
  m_mesh_fields[ NODE_CENTERED ]->resize( getNumberOfNodes() );
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
  return m_mesh_fields[ NODE_CENTERED ]->checkConsistency( getNumberOfNodes(),
                                                           getNodeCapacity() );
}

} /* namespace mint */
} /* namespace axom */
