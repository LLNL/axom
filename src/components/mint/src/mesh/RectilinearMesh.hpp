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

#ifndef RECTILINEARMESH_HXX_
#define RECTILINEARMESH_HXX_

#include "axom/Macros.hpp"
#include "axom/Types.hpp"

#include "mint/Array.hpp"
#include "mint/MeshCoordinates.hpp"
#include "mint/StructuredMesh.hpp"
#include "mint/config.hpp"

namespace axom
{
namespace mint
{

class RectilinearMesh : public StructuredMesh
{
public:


  RectilinearMesh() = delete;

  /*!
   * \brief Constructs a rectilinear mesh instance.
   * \param [in] dimension the dimension of the mesh.
   * \param [in] ext the mesh extent.
   */
  RectilinearMesh( int dimension, int64 ext[6] );

  /*!
   * \brief Destructor.
   */
  virtual ~RectilinearMesh();

  /*!
   * \brief Sets the coordinate along the given dimension.
   * \param [in] dim the dimension in query.
   * \param [in] i the index of the coordinate to set.
   * \param [in] coord the coordinate to set.
   * \pre dim >= 0 && dim < this->getDimension()
   * \pre i >= 0 && i < ndims[ i ]
   */
  inline void setCoordinate( int dim, IndexType i, double coord ) const;

  /*!
   * \brief Gets the coordinate along the given dimension.
   * \return the position of the ith coordinate along the given dimension.
   * \pre dim >= 0 && dim < this->getDimension()
   * \pre i >= 0 && i < ndims[ i ]
   */
  inline double getCoordinate( int dim, IndexType i ) const;


  virtual double* getCoordinateArray( int dim ) final override;

  virtual const double* getCoordinateArray( int dim ) const final override;

  /*!
   * \brief Copy the coordinates of the given node into the provided buffer.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] coords the buffer to copy the coordinates into, of length at
   *  least getDimension().
   *
   * \note provided only for convenience, do not use inside a loop. Instead use
   *  getCoordinate() or getCoordinateArray() methods to calculate the nodal 
   *  coordinates.
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != AXOM_NULLPTR
   */
  virtual void getNode( IndexType nodeID, double* node ) const override final;

private:

  Array< double >* m_coordinates[3];

  DISABLE_COPY_AND_ASSIGNMENT(RectilinearMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(RectilinearMesh);
};

//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void RectilinearMesh::setCoordinate( int dim, IndexType i,
                                            double coord ) const
{
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );
  SLIC_ASSERT( i >= 0 && i < m_coordinates[ dim ]->size() );

  (*m_coordinates[ dim ])( i ) = coord;
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getCoordinate( int dim, IndexType i ) const
{
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );
  SLIC_ASSERT( i >= 0 && i < m_coordinates[ dim ]->size() );

  return (*m_coordinates[ dim ])( i );
}

//------------------------------------------------------------------------------
inline double* RectilinearMesh::getCoordinateArray( int dim )
{
  SLIC_ASSERT( 0 <= dim && dim < getDimension() );
  SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
inline const double* RectilinearMesh::getCoordinateArray( int dim ) const
{
  SLIC_ASSERT( 0 <= dim && dim < getDimension() );
  SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
inline void RectilinearMesh::getNode( IndexType nodeID, double* node ) const
{
  SLIC_ASSERT( 0 <= nodeID && nodeID < getNumberOfNodes() );
  SLIC_ASSERT( node != AXOM_NULLPTR );

  const int n_dims = getDimension();
  IndexType gridIndices[3];

  if ( n_dims == 1 )
  {
    node[ 0 ] = getCoordinate( 0, nodeID );
    return;
  }
  else if ( n_dims == 2 )
  {
    m_extent.getGridIndex( nodeID, gridIndices[0], gridIndices[1] );
  }
  else
  {
    SLIC_ASSERT( n_dims == 3 );
    m_extent.getGridIndex( nodeID, gridIndices[0], gridIndices[1], gridIndices[2] );
  }

  for ( int dim = 0; dim < n_dims; ++dim )
  {
    node[ dim ] = getCoordinate( dim, gridIndices[ dim ] );
  }
}

} /* namespace mint */
} /* namespace axom */

#endif /* RECTILINEARMESH_HXX_ */
