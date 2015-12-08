/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file ParticleMesh.hxx
 *
 * \date Sep 27, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */


#include "quest/ParticleMesh.hpp"
#include "common/CommonTypes.hpp"

namespace meshtk
{

ParticleMesh::ParticleMesh( ) :
    Mesh( -1, UNDEFINED_MESH, -1, -1 ),
    m_particle_coordinates( ATK_NULLPTR )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension ) :
    Mesh( dimension, PARTICLE_MESH, 0, 0 ),
    m_particle_coordinates( new MeshCoordinates(dimension) )
{

}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, int blockId, int partId ) :
    Mesh( dimension, PARTICLE_MESH, blockId, partId ),
    m_particle_coordinates( new MeshCoordinates(dimension) )
{

}

//------------------------------------------------------------------------------
ParticleMesh::~ParticleMesh()
{
  delete m_particle_coordinates;
  m_particle_coordinates = ATK_NULLPTR;
}

} /* namespace meshtk */
