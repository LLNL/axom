/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


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


#include "mint/ParticleMesh.hpp"
#include "common/CommonTypes.hpp"

namespace mint
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

} /* namespace mint */
