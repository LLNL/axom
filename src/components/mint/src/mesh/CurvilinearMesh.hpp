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

#ifndef CURVILINEARMESH_HXX_
#define CURVILINEARMESH_HXX_

#include "mint/StructuredMesh.hpp"
#include "mint/MeshCoordinates.hpp"
#include "slic/slic.hpp"

namespace axom
{
namespace mint
{

class CurvilinearMesh : public StructuredMesh
{
public:

  /*!
   * \brief Constructs a curvilinear mesh instance.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the logical extent of this mesh instance.
   */
  CurvilinearMesh( int dimension, const globalIndex ext[6] );

  /*!
   * \brief Constructs a curvilinear mesh instance.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the logical extent of this mesh instance.
   * \param [in] blockId the block ID of this mesh
   * \param [in] partId the partition ID of this mesh
   */
  CurvilinearMesh( int dimension, const globalIndex ext[6], int blockId,
                   int partId );

  /*!
   * \brief Destructor.
   */
  virtual ~CurvilinearMesh()
  {}

  /// \name GetNode() methods
  /// @{

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre coordinates != AXOM_NULLPTR.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes().
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
   * \brief Returns the coordinates of the node at (i,j)
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
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes()
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

  /// \name SetNode() methods
  /// @{

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \param [in] z the z--coorindate to set at the given node.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes
   * \pre this->getDimension() == 3
   */
  void setNode( localIndex nodeIdx, double x, double y, double z );

  /*!
   * \brief Sets the coordinates of the node at (i,j,k).
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \param [in] z the z--coordinate to set at the given node.
   * \pre this->getDimension() == 3
   */
  void setNode( localIndex i, localIndex j, localIndex k, double x, double y,
                double z );

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \pre this->getDimension() == 2
   */
  void setNode( localIndex nodeIdx, double x, double y );

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \pre this->getDimension() == 2
   */
  void setNode( localIndex i, localIndex j, double x, double y );

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \pre this->getDimension() == 1
   */
  void setNode( localIndex nodeIdx, double x );

  /// @}

  /*!
   * \brief Returns a pointer to the coordinates array for the given dimension.
   * \param [in] dim the requested coordinate dimension.
   * \return ptr pointer to the coordinates array.
   * \pre dim >= 0 && dim < this->getDimension().
   * \post ptr != AXOM_NULLPTR.
   */
  const double* getMeshCoordinateArray( int dim ) const;

private:

  MeshCoordinates m_coordinates;

  CurvilinearMesh();

  CurvilinearMesh(const CurvilinearMesh&); // Not implemented
  CurvilinearMesh& operator=(const CurvilinearMesh&); // Not implemented
};


//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline const double* CurvilinearMesh::getMeshCoordinateArray( int dim ) const
{
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );

  return( m_coordinates.getCoordinateArray( dim ) );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( localIndex nodeIdx, double x, double y,
                                      double z)
{
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );
  SLIC_ASSERT( this->getDimension() == 3 );

  m_coordinates.setPoint( nodeIdx, x, y, z );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( localIndex i, localIndex j, localIndex k,
                                      double x, double y, double z )
{
  SLIC_ASSERT( this->getDimension() == 3 );

  localIndex nodeIdx = m_extent.getLinearIndex( i, j, k );
  this->setNode( nodeIdx, x, y, z );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( localIndex nodeIdx, double x, double y )
{
  SLIC_ASSERT(  nodeIdx >= 0 && nodeIdx < this->getMeshNumberOfNodes() );
  SLIC_ASSERT(  this->getDimension() == 2 );

  m_coordinates.setPoint( nodeIdx, x, y );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( localIndex i, localIndex j, double x,
                                      double y )
{
  SLIC_ASSERT( this->getDimension() == 2 );

  localIndex nodeIdx = m_extent.getLinearIndex( i, j );
  this->setNode( nodeIdx, x, y );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( localIndex nodeIdx, double x )
{
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < this->getMeshNumberOfNodes() );
  SLIC_ASSERT( this->getDimension() == 1 );

  m_coordinates.setPoint( nodeIdx, x );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::getNode( localIndex nodeIdx,
                                      double* coordinates) const
{
  SLIC_ASSERT( coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );

  for ( int i=0 ; i < this->getDimension() ; ++i )
  {
    const double* coords = m_coordinates.getCoordinateArray( i );
    coordinates[ i ] = coords[ nodeIdx ];
  }
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::getNode( localIndex i, localIndex j,
                                      double* coordinates ) const
{
  SLIC_ASSERT( coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( this->getDimension() == 2 );

  const localIndex nodeIdx = m_extent.getLinearIndex( i, j );
  this->getNode( nodeIdx, coordinates );
}

//------------------------------------------------------------------------------
inline
void CurvilinearMesh::getNode( localIndex i, localIndex j, localIndex k,
                               double* coordinates) const
{
  SLIC_ASSERT( coordinates != AXOM_NULLPTR );
  SLIC_ASSERT( this->getDimension() == 3 );

  const localIndex nodeIdx = m_extent.getLinearIndex( i, j, k );
  this->getNode( nodeIdx, coordinates );
}

//------------------------------------------------------------------------------
inline double CurvilinearMesh::getNodeCoordinate( localIndex nodeIdx,
                                                  int dim ) const
{
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );
  SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );

  const double* coord = m_coordinates.getCoordinateArray( dim );
  return( coord[nodeIdx] );
}

//------------------------------------------------------------------------------
inline double CurvilinearMesh::getNodeCoordinate( localIndex i, localIndex j,
                                                  int dim ) const
{
  SLIC_ASSERT( this->getDimension() == 2 );

  const localIndex nodeIdx = m_extent.getLinearIndex( i, j );
  return ( this->getNodeCoordinate(nodeIdx, dim) );
}

//------------------------------------------------------------------------------
inline
double CurvilinearMesh::getNodeCoordinate( localIndex i, localIndex j,
                                           localIndex k, int dim ) const
{
  SLIC_ASSERT( this->getDimension() == 3 );

  const localIndex nodeIdx = m_extent.getLinearIndex( i, j, k );
  return ( this->getNodeCoordinate(nodeIdx, dim) );
}

} /* namespace mint */
} /* namespace axom */

#endif /* CURVILINEARMESH_HXX_ */
