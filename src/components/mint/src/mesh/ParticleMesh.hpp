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

#ifndef PARTICLEMESH_HXX_
#define PARTICLEMESH_HXX_

#include "mint/Mesh.hpp"
#include "mint/MeshCoordinates.hpp"
#include "slic/slic.hpp"
#include "mint/config.hpp"           /* For IndexType */

namespace axom
{
namespace mint
{

class ParticleMesh : public Mesh
{

public:

  /*!
   * \brief Constructs a ParticleMesh instance.
   * \param [in] dimension the ambient dimension of the particle mesh.
   */
  ParticleMesh( int dimension, IndexType particleCapacity );

  /*!
   * \brief Constructs a ParticleMesh instance.
   * \param [in] dimension the ambient dimension of the particle mesh.
   * \param [in] blockId the block ID.
   * \param [in] partId the partition ID.
   */
  ParticleMesh( int dimension, IndexType particleCapacity, int blockId,
                int partId );

  /*!
   * \brief Destructor.
   */
  virtual ~ParticleMesh() { }

  /*!
   * \brief Returns the total number of particles in this mesh instance.
   * \return N the total number of particles.
   */
  IndexType getNumberOfParticles() const
  { return m_particle_coordinates.numNodes(); }

  /*!
   * \brief Returns the particle coordinates array for the given dimension.
   * \param [in] idim the requested dimension.
   * \return coords pointer to the coorindates array.
   * \pre idim >= 0 && idim < this->getDimension()
   */
  double* getParticlesCoordinatesArray( int idim );

  /*!
   * \brief Returns the particle coordinates array for the given dimension.
   * \param [in] idim the requested dimension.
   * \return coords pointer to the coorindates array.
   * \pre idim >= 0 && idim < this->getDimension()
   */
  const double* getParticlesCoordinatesArray( int idim ) const;

  /*!
   * \brief Add a particle in to the particle mesh.
   * \param [in] x the x--coordinate of the particle.
   * \pre this->getDimension() == 1
   */
  void addParticle( double x );

  /*!
   * \brief Add a particle in to the particle mesh.
   * \param [in] x the x--coordinate of the particle.
   * \param [in] y the y--coordinate of the particle.
   * \pre this->getDimension() == 2
   */
  void addParticle( double x, double y );

  /*!
   * \brief Add a particle in this particle mesh instance.
   * \param [in] x the x--coordinate of the particle.
   * \param [in] y the y--coordinate of the particle.
   * \param [in] z the z--coordinate of the particle.
   * \pre this->getDimension() == 3.
   */
  void addParticle( double x, double y, double z );

  /*!
   * \brief Gets the particle particle coordinates of a given particle.
   * \param [in] partIdx the index of the particle in question.
   * \param [in] part_coords user-supplied buffer for the particle coordinates.
   * \pre partIdx >= 0 && partIdx < this->getNumberOfParticles()
   */
  void getParticleCoordinates( IndexType partIdx,
                               double part_coords[3] ) const;


  IndexType getCapacity() const
  { return m_particle_coordinates.capacity(); }


  void reserve( IndexType capacity )
  { return m_particle_coordinates.reserve( capacity ); }


  double getResizeRatio() const
  { return m_particle_coordinates.getResizeRatio(); }


private:

  MeshCoordinates m_particle_coordinates;

  DISABLE_COPY_AND_ASSIGNMENT(ParticleMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(ParticleMesh);
};


//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline double* ParticleMesh::getParticlesCoordinatesArray( int idim )
{
  SLIC_ASSERT( idim >= 0 && idim < this->getDimension() );
  return m_particle_coordinates.getCoordinateArray( idim );
}

//------------------------------------------------------------------------------
inline const double* ParticleMesh::getParticlesCoordinatesArray( int idim )
const
{
  SLIC_ASSERT( idim >= 0 && idim < this->getDimension() );
  return m_particle_coordinates.getCoordinateArray( idim );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::addParticle( double x ) {
  SLIC_ASSERT( this->getDimension() == 1 );
  m_particle_coordinates.append( &x );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::addParticle( double x, double y ) {
  SLIC_ASSERT( this->getDimension() == 2 );
  double xx[2] = { x, y };
  m_particle_coordinates.append( xx );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::addParticle( double x, double y, double z ) {
  SLIC_ASSERT( this->getDimension() == 3);
  double xx[3] = { x, y, z };
  m_particle_coordinates.append( xx );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::getParticleCoordinates( IndexType partIdx,
                                                  double part_coords[3] ) const
{
  SLIC_ASSERT( partIdx >= 0 && partIdx < this->getNumberOfParticles() );

  for ( int i=0 ; i < this->getDimension() ; ++i )
  {
    const double* px = this->getParticlesCoordinatesArray( i );
    part_coords[ i ] = px[ partIdx ];
  }

}

} /* namespace mint */
} /* namespace axom */

#endif /* PARTICLEMESH_HXX_ */
