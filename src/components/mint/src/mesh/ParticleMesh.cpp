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
#include "axom/Types.hpp"

namespace axom
{
namespace mint
{

ParticleMesh::ParticleMesh( ) :
  Mesh( -1, MINT_UNDEFINED_MESH, -1, -1 ),
  m_particle_coordinates( AXOM_NULLPTR )
{}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, int particleCapacity ):
  Mesh( dimension, MINT_PARTICLE_MESH, 0, 0 ),
  m_particle_coordinates( new MeshCoordinates( dimension, particleCapacity) )
{}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, int particleCapacity, int blockId, 
																 int partId ):
  Mesh( dimension, MINT_PARTICLE_MESH, blockId, partId ),
  m_particle_coordinates( new MeshCoordinates( dimension, particleCapacity ) )
{}

//------------------------------------------------------------------------------
ParticleMesh::~ParticleMesh()
{
  delete m_particle_coordinates;
  m_particle_coordinates = AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
