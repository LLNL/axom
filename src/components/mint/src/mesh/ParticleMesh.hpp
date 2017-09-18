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
#include "mint/CellType.hpp"
#include "mint/MeshCoordinates.hpp"
#include "slic/slic.hpp"

#include <cstddef> // for AXOM_NULLPTR

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
  explicit ParticleMesh( int dimension, int particleCapacity );

  /*!
   * \brief Constructs a ParticleMesh instance.
   * \param [in] dimension the ambient dimension of the particle mesh.
   * \param [in] blockId the block ID.
   * \param [in] partId the partition ID.
   */
  ParticleMesh( int dimension, int particleCapacity, int blockId, int partId );

  /*!
   * \brief Destructor.
   */
  virtual ~ParticleMesh();

  /// \name Virtual Mesh API
  /// @{

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfNodes() const
  { return this->getNumberOfParticles(); };

  /*!
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfCells() const
  { return this->getNumberOfParticles(); };

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual int getMeshNumberOfCellNodes( int /*cellIdx*/ ) const { return 1; };

  /*!
   * \brief Returns the cell connectivity of the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \param [out] cell user-supplied buffer to store cell connectivity info.
   * \note cell must have sufficient size to hold the connectivity information.
   * \pre cellIdx >= 0 && cellIdx < this->getMeshNumberOfCells()
   * \pre cell != AXOM_NULLPTR.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshCell( int cellIdx, int* cell ) const {
    cell[0]=cellIdx;
  };

  /*!
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   */
  virtual int getMeshCellType( int /*cellIdx*/ ) const { return MINT_VERTEX; };

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] coordinates user-supplied buffer to store the node coordinates.
   * \pre nodeIdx >= && nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshNode( int nodeIdx, double* coordinates ) const
  { this->getParticleCoordinates( nodeIdx, coordinates ); };

  /*!
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < m_ndims
   */
  virtual double getMeshNodeCoordinate( int nodeIdx, int dim ) const
  { return this->getParticlesCoordinatesArray(dim)[ nodeIdx ]; };

  /// @}

  /*!
   * \brief Returns the total number of particles in this mesh instance.
   * \return N the total number of particles.
   */
  int getNumberOfParticles() const
  { return m_particle_coordinates->getSize(); }

  /*!
   * \brief Returns the particle coordinates array for the given dimension.
   * \param [in] idim the requested dimension.
   * \return coords pointer to the coorindates array.
   * \pre idim >= 0 && idim < this->getDimension()
   */
  double* getParticlesCoordinatesArray( int idim ) const;

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
  void getParticleCoordinates( int partIdx, double part_coords[3] ) const;

private:

  /*!
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent users from calling it.
   */
  ParticleMesh();

  MeshCoordinates* m_particle_coordinates;

  DISABLE_COPY_AND_ASSIGNMENT(ParticleMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(ParticleMesh);
};

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------
namespace axom
{
namespace mint
{

inline double* ParticleMesh::getParticlesCoordinatesArray( int idim ) const
{
  SLIC_ASSERT( idim >= 0 && idim < this->getDimension() );
  return m_particle_coordinates->getCoordinateArray( idim );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::addParticle( double x ) {
  SLIC_ASSERT( this->getDimension() == 1 );
  m_particle_coordinates->addPoint( x );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::addParticle( double x, double y ) {
  SLIC_ASSERT( this->getDimension() == 2 );
  m_particle_coordinates->addPoint( x, y );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::addParticle( double x, double y, double z ) {
  SLIC_ASSERT( this->getDimension() == 3);
  m_particle_coordinates->addPoint( x, y, z );
}

//------------------------------------------------------------------------------
inline void ParticleMesh::getParticleCoordinates( int partIdx,
                                                  double part_coords[3] ) const
{
  SLIC_ASSERT( partIdx >= 0 && partIdx < this->getNumberOfParticles() );

  for ( int i=0 ; i < this->getDimension() ; ++i )
  {
    double* px = this->getParticlesCoordinatesArray( i );
    part_coords[ i ] = px[ partIdx ];
  }

}

} /* namespace mint */
} /* namespace axom */

#endif /* PARTICLEMESH_HXX_ */
