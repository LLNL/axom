/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#include "mint/ParticleMesh.hpp"
#include "common/CommonTypes.hpp"

namespace axom {
namespace mint {

ParticleMesh::ParticleMesh( ) :
    Mesh( -1, MINT_UNDEFINED_MESH, -1, -1 ),
    m_particle_coordinates( ATK_NULLPTR )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension ) :
    Mesh( dimension, MINT_PARTICLE_MESH, 0, 0 ),
    m_particle_coordinates( new MeshCoordinates(dimension) )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, int blockId, int partId ) :
    Mesh( dimension, MINT_PARTICLE_MESH, blockId, partId ),
    m_particle_coordinates( new MeshCoordinates(dimension) )
{

}

//------------------------------------------------------------------------------
ParticleMesh::~ParticleMesh()
{
  delete m_particle_coordinates;
  m_particle_coordinates = ATK_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */
