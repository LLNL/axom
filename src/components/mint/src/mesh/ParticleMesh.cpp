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
#include "mint/config.hpp"
#include "axom/Types.hpp"

namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, IndexType particleCapacity ) :
  Mesh( dimension, MINT_PARTICLE_MESH, 0, 0 ),
  m_particle_coordinates( dimension, particleCapacity, particleCapacity, 2.0 )
{}

//------------------------------------------------------------------------------
ParticleMesh::ParticleMesh( int dimension, IndexType particleCapacity,
                            int blockId, int partId ) :
  Mesh( dimension, MINT_PARTICLE_MESH, blockId, partId ),
  m_particle_coordinates( dimension, particleCapacity, particleCapacity, 2.0 )
{}

} /* namespace mint */
} /* namespace axom */
