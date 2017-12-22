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

#include "mint/StructuredMesh.hpp"
#include "mint/MeshCoordinates.hpp"
#include "mint/Array.hpp"
#include "mint/DataTypes.hpp"

#include "axom/Macros.hpp"
#include "axom/Types.hpp"

#define X_COORDINATE 0
#define Y_COORDINATE 1
#define Z_COORDINATE 2

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
  RectilinearMesh( int dimension, globalIndex ext[6] );

  /*!
   * \brief Constructs a rectilinear mesh instance.
   * \param [in] dimension the dimension of the mesh.
   * \param [in] ext the mesh extent.
   * \param [in] blockId the block ID.
   * \param [in] partitionId the partition ID.
   */
  RectilinearMesh( int dimension, globalIndex ext[6], int blockId,
                   int partitionId );

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
  inline void setCoordinate( int dim, localIndex i, double coord ) const;

  /*!
   * \brief Gets the coordinate along the given dimension.
   * \return the position of the ith coordinate along the given dimension.
   * \pre dim >= 0 && dim < this->getDimension()
   * \pre i >= 0 && i < ndims[ i ]
   */
  inline double getCoordinate( int dim, localIndex i ) const;

  /*!
   * \brief Returns pointer to the coordinates array along the given dimension.
   * \param [in] dim the requested dimension.
   * \return coordsPtr pointer to the coordinates array.
   */
  const double* getCoordinateArray( int dim ) const;

  /// \name GetNode() methods.
  /// @{

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre coordinates != AXOM_NULLPTR.
   * \pre nodeIdx >= 0 && nodeIdx < getNumberOfNodes().
   */
  virtual void getNode( localIndex nodeIdx,
                        double* coordinates ) const override;

  /*!
   * \brief Returns the coordinates of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 2
   */
  virtual void getNode( localIndex i, localIndex j,
                        double* coordinates ) const override;

  /*!
   * \brief Returns the coordinates of the node at (i,j,k)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 3
   */
  virtual void getNode( localIndex i, localIndex j, localIndex k,
                        double* coordinates ) const override;

  /*!
   * \brief Returns the coordinate of the given node.
   * \param [in] nodeIdx index of the node in query.
   * \param [in] dim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre nodeIdx >= 0 && nodeIdx < getNumberOfNodes()
   * \pre dim >= 0 && dim < m_ndims.
   */
  virtual double getNodeCoordinate( localIndex nodeIdx,
                                    int dim ) const override;

  /*!
   * \brief Returns the coordinate value of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] dim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==2.
   * \pre dim >= 0 && dim < m_ndims.
   */
  virtual double getNodeCoordinate( localIndex i, localIndex j,
                                    int dim ) const override;

  /*!
   * \brief Returns the coordinate value of the node at (i,j,k)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [in] dim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==3.
   * \pre dim >= 0 && dim < m_ndims.
   */
  virtual double getNodeCoordinate( localIndex i, localIndex j, localIndex k,
                                    int dim ) const override;

  /// @}

private:

  Array< double >* m_coordinates[3];

  DISABLE_COPY_AND_ASSIGNMENT(RectilinearMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(RectilinearMesh);
};

//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void RectilinearMesh::setCoordinate( int dim, localIndex i,
                                            double coord ) const
{
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );
  SLIC_ASSERT( i >= 0 && i < m_coordinates[ dim ]->size() );

  (*m_coordinates[ dim ])( i ) = coord;
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getCoordinate( int dim, localIndex i ) const
{
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );
  SLIC_ASSERT( i >= 0 && i < m_coordinates[ dim ]->size() );

  return (*m_coordinates[ dim ])( i );
}

//------------------------------------------------------------------------------
inline const double* RectilinearMesh::getCoordinateArray( int dim ) const
{
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
inline void RectilinearMesh::getNode( localIndex nodeIdx,
                                      double* coordinates ) const
{
  SLIC_ASSERT( coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < getNumberOfNodes() );

  localIndex ijk[3];
  m_extent.getGridIndex( nodeIdx, ijk[0], ijk[1], ijk[2] );

  for ( int dim = 0 ; dim < this->getDimension() ; ++dim )
  { coordinates[ dim ] = (*m_coordinates[ dim ])( ijk[ dim ] ); }
}

//------------------------------------------------------------------------------
inline void RectilinearMesh::getNode( localIndex i, localIndex j,
                                      double* coordinates ) const
{
  SLIC_ASSERT( coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( this->getDimension() == 2 );
  SLIC_ASSERT( i >= 0 && j >= 0 );

  coordinates[ 0 ] = (*m_coordinates[ X_COORDINATE ])( i );
  coordinates[ 1 ] = (*m_coordinates[ Y_COORDINATE ])( j );
}

//------------------------------------------------------------------------------
inline void RectilinearMesh::getNode( localIndex i, localIndex j, localIndex k,
                                      double* coordinates ) const
{
  SLIC_ASSERT( coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( this->getDimension() == 3 );
  SLIC_ASSERT( i >= 0 && j >= 0 && k >= 0 );

  coordinates[ 0 ] = (*m_coordinates[ X_COORDINATE ])( i );
  coordinates[ 1 ] = (*m_coordinates[ Y_COORDINATE ])( j );
  coordinates[ 2 ] = (*m_coordinates[ Z_COORDINATE ])( k );
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getNodeCoordinate( localIndex nodeIdx,
                                                  int dim ) const
{
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < getNumberOfNodes() );
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );

  localIndex ijk[3];
  m_extent.getGridIndex( nodeIdx, ijk[0], ijk[1], ijk[2] );
  return (*m_coordinates[ dim ])( ijk[dim]  );
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getNodeCoordinate( localIndex i, localIndex j,
                                                  int dim ) const
{
  SLIC_ASSERT( this->getDimension() == 2 );
  SLIC_ASSERT( dim >= 0 && dim < 2 );

  localIndex ijk[2] = { i, j };
  return (*m_coordinates[ dim ])( ijk[dim] );
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getNodeCoordinate( localIndex i, localIndex j,
                                                  localIndex k, int dim ) const
{
  SLIC_ASSERT( this->getDimension()==3 );
  SLIC_ASSERT( dim >= 0 && dim < 3 );

  localIndex ijk[3] = { i, j, k };
  return (*m_coordinates[ dim ])( ijk[dim] );
}

} /* namespace mint */
} /* namespace axom */

#endif /* RECTILINEARMESH_HXX_ */
