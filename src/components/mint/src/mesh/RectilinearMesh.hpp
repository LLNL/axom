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

#include "axom/Macros.hpp"
#include "axom/Types.hpp"

namespace axom
{
namespace mint
{

class RectilinearMesh : public StructuredMesh
{
public:

  /*!
   * \brief Constructs a rectilinear mesh instance.
   * \param [in] dimension the dimension of the mesh.
   * \param [in] ext the mesh extent.
   */
  RectilinearMesh( int dimension, int ext[6] );

  /*!
   * \brief Constructs a rectilinear mesh instance.
   * \param [in] dimension the dimension of the mesh.
   * \param [in] ext the mesh extent.
   * \param [in] blockId the block ID.
   * \param [in] partitionId the partition ID.
   */
  RectilinearMesh( int dimension, int ext[6], int blockId, int partitionId );

  /*!
   * \brief Destructor.
   */
  virtual ~RectilinearMesh();

  /*!
   * \brief Sets the coordinate along the given dimension.
   * \param [in] idim the dimension in query.
   * \param [in] i the index of the coordinate to set.
   * \param [in] coord the coordinate to set.
   * \pre idim >= 0 && idim < this->getDimension()
   * \pre i >= 0 && i < ndims[ i ]
   */
  void setCoordinate( int idim, int i, double coord );

  /*!
   * \brief Returns pointer to the coordinates array along the given dimension.
   * \param [in] idim the requested dimension.
   * \return coordsPtr pointer to the coordinates array.
   */
  const double* getCoordinateArray( int idim ) const
  { return m_coordinates->getCoordinateArray( idim ); };

  /// \name GetNode() methods.
  /// @{

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre coordinates != AXOM_NULLPTR.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes().
   */
  virtual void getNode( int nodeIdx, double* coordinates ) const;

  /*!
   * \brief Returns the coordinates of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 2
   */
  virtual void getNode( int i, int j, double* coordinates ) const;

  /*!
   * \brief Returns the coordinates of the node at (i,j,k)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 3
   */
  virtual void getNode( int i, int j, int k, double* coordinates ) const;

  /*!
   * \brief Returns the coordinate of the given node.
   * \param [in] nodeIdx index of the node in query.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes()
   * \pre idim >= 0 && idim < m_ndims.
   */
  virtual double getNodeCoordinate( int nodeIdx, int idim  ) const;

  /*!
   * \brief Returns the coordinate value of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==2.
   * \pre idim >= 0 && idim < m_ndims.
   */
  virtual double getNodeCoordinate( int i, int j, int idim ) const;

  /*!
   * \brief Returns the coordinate value of the node at (i,j,k)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==3.
   * \pre idim >= 0 && idim < m_ndims.
   */
  virtual double getNodeCoordinate( int i, int j, int k, int idim ) const;

  /// @}

private:

  /*!
   * \brief Default constructor.
   * \note Made private to prevent users from calling it.
   */
  RectilinearMesh();

  MeshCoordinates* m_coordinates;

  DISABLE_COPY_AND_ASSIGNMENT(RectilinearMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(RectilinearMesh);
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

inline void RectilinearMesh::setCoordinate( int idim, int i, double coord )
{
  SLIC_ASSERT(  idim >= 0 && idim < this->getDimension() );
  SLIC_ASSERT(  i >= 0 && i < m_coordinates->getCoordinateArraySize( idim ) );

  double* xc = m_coordinates->getCoordinateArray( idim );
  xc[ i ]    = coord;
}

//------------------------------------------------------------------------------
inline void RectilinearMesh::getNode( int nodeIdx, double* coordinates ) const
{
  SLIC_ASSERT(  coordinates != AXOM_NULLPTR );
  SLIC_ASSERT(  nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );

  int ijk[3];
  m_extent->getGridIndex( nodeIdx, ijk[0], ijk[1], ijk[2] );

  for ( int i=0 ; i < this->getDimension() ; ++i )
  {
    const double* xc = this->getCoordinateArray( i );
    coordinates[ i ] = xc[ ijk[i] ];
  }

}

//------------------------------------------------------------------------------
inline void RectilinearMesh::getNode( int i, int j, double* coordinates ) const
{
  SLIC_ASSERT(  coordinates != AXOM_NULLPTR );
  SLIC_ASSERT(  this->getDimension()==2 );

  const double* xc = this->getCoordinateArray( 0 );
  const double* yc = this->getCoordinateArray( 1 );
  coordinates[ 0 ] = xc[ i ];
  coordinates[ 1 ] = yc[ j ];
}

//------------------------------------------------------------------------------
inline void RectilinearMesh::getNode(
  int i, int j, int k, double* coordinates ) const
{
  SLIC_ASSERT(  coordinates != AXOM_NULLPTR );
  SLIC_ASSERT(  this->getDimension()==3 );

  const double* xc = this->getCoordinateArray( 0 );
  const double* yc = this->getCoordinateArray( 1 );
  const double* zc = this->getCoordinateArray( 2 );
  coordinates[ 0 ] = xc[ i ];
  coordinates[ 1 ] = yc[ j ];
  coordinates[ 2 ] = zc[ k ];
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getNodeCoordinate( int nodeIdx, int idim  ) const
{
  SLIC_ASSERT(  nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );
  SLIC_ASSERT(  idim >= 0 && idim < this->getDimension() );

  int ijk[3];
  m_extent->getGridIndex( nodeIdx, ijk[0], ijk[1], ijk[2] );

  const double* xc = this->getCoordinateArray( idim );
  return xc[ ijk[idim]  ];
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getNodeCoordinate( int i, int j, int idim ) const
{
  SLIC_ASSERT(  this->getDimension()==2 );
  SLIC_ASSERT(  idim >= 0 && idim < 2 );

  int ijk[2] = { i, j };

  const double* xc = this->getCoordinateArray( idim );
  return xc[ ijk[idim] ];
}

//------------------------------------------------------------------------------
inline double RectilinearMesh::getNodeCoordinate(
  int i, int j, int k, int idim ) const
{
  SLIC_ASSERT(  this->getDimension()==3 );
  SLIC_ASSERT(  idim >= 0 && idim < 3 );

  int ijk[3] = { i, j, k };

  const double* xc = this->getCoordinateArray( idim );
  return xc[ ijk[idim] ];
}

} /* namespace mint */
} /* namespace axom */
#endif /* RECTILINEARMESH_HXX_ */
