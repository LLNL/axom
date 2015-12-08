/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file CurvilinearMesh.hxx
 *
 * \date Sep 20, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef CURVILINEARMESH_HXX_
#define CURVILINEARMESH_HXX_

#include "quest/StructuredMesh.hpp"
#include "quest/MeshCoordinates.hpp"

namespace meshtk
{

class CurvilinearMesh : public StructuredMesh
{
public:

  /*!
   *****************************************************************************
   * \brief Constructs a curvilinear mesh instance.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the logical extent of this mesh instance.
   *****************************************************************************
   */
  CurvilinearMesh( int dimension, int ext[6] );

  /*!
   *****************************************************************************
   * \brief Constructs a curvilinear mesh instance.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the logical extent of this mesh instance.
   * \param [in] blockId the block ID of this mesh
   * \param [in] partId the partition ID of this mesh
   *****************************************************************************
   */
  CurvilinearMesh( int dimension, int ext[6], int blockId, int partId );

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  virtual ~CurvilinearMesh();

  /// \name GetNode() methods
  /// @{

  /*!
   *****************************************************************************
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre coordinates != ATK_NULLPTR.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes().
   *****************************************************************************
   */
  virtual void getNode( int nodeIdx, double* coordinates ) const;

  /*!
   *****************************************************************************
   * \brief Returns the coordinates of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 2
   *****************************************************************************
   */
  virtual void getNode( int i, int j, double* coordinates ) const;

  /*!
   *****************************************************************************
   * \brief Returns the coordinates of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 3
   *****************************************************************************
   */
  virtual void getNode( int i, int j, int k, double* coordinates ) const;

  /*!
   *****************************************************************************
   * \brief Returns the coordinate of the given node.
   * \param [in] nodeIdx index of the node in query.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes()
   * \pre idim >= 0 && idim < m_ndims.
   *****************************************************************************
   */
  virtual double getNodeCoordinate( int nodeIdx, int idim  ) const;

  /*!
   *****************************************************************************
   * \brief Returns the coordinate value of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==2.
   * \pre idim >= 0 && idim < m_ndims.
   *****************************************************************************
   */
  virtual double getNodeCoordinate( int i, int j, int idim ) const;

  /*!
   *****************************************************************************
   * \brief Returns the coordinate value of the node at (i,j,k)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==3.
   * \pre idim >= 0 && idim < m_ndims.
   *****************************************************************************
   */
  virtual double getNodeCoordinate( int i, int j, int k, int idim ) const;

  /// @}

  /// \name SetNode() methods
  /// @{

  /*!
   *****************************************************************************
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \param [in] z the z--coorindate to set at the given node.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes
   * \pre this->getDimension()==3
   *****************************************************************************
   */
  void setNode( int nodeIdx, double x, double y, double z );

  /*!
   *****************************************************************************
   * \brief Sets the coordinates of the node at (i,j,k).
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \param [in] z the z--coordinate to set at the given node.
   * \pre this->getDimension()==3
   *****************************************************************************
   */
  void setNode( int i, int j, int k,
                double x, double y, double z );

  /*!
   *****************************************************************************
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \pre this->getDimension()==2
   *****************************************************************************
   */
  void setNode( int nodeIdx, double x, double y );


  /*!
   *****************************************************************************
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \pre this->getDimension()==2
   *****************************************************************************
   */
  void setNode( int i, int j, double x, double y );

  /// @}

  /*!
   *****************************************************************************
   * \brief Returns a pointer to the coordinates array for the given dimension.
   * \param [in] idim the requested coordinate dimension.
   * \return ptr pointer to the coordinates array.
   * \pre idim >= 0 && idim < this->getDimension().
   * \post ptr != ATK_NULLPTR.
   *****************************************************************************
   */
  const double* getMeshCoordinateArray( int idim ) const;

private:

  /*!
   ***************************************************************************
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent applications from calling it.
   ***************************************************************************
   */
  CurvilinearMesh();

  MeshCoordinates* m_coordinates;

  CurvilinearMesh(const CurvilinearMesh&); // Not implemented
  CurvilinearMesh& operator=(const CurvilinearMesh&); // Not implemented
};

} /* namespace meshtk */

//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

namespace meshtk
{


inline const double* CurvilinearMesh::getMeshCoordinateArray( int idim ) const
{
  assert( idim >= 0 && idim < this->getDimension() );
  return( m_coordinates->getCoordinateArray( idim ) );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( int nodeIdx, double x, double y, double z)
{
  assert( nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );
  assert( this->getDimension()==3 );
  m_coordinates->setPoint( nodeIdx, x, y, z );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( int i,    int j,    int k,
                                     double x, double y, double z )
{
  int nodeIdx = m_extent->getLinearIndex( i, j, k );
  this->setNode( nodeIdx, x, y, z );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( int nodeIdx, double x, double y )
{
  assert( nodeIdx >= 0 && nodeIdx < this->getMeshNumberOfNodes() );
  assert( this->getDimension()==2 );
  m_coordinates->setPoint( nodeIdx, x, y );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( int i, int j, double x, double y )
{
  int nodeIdx = m_extent->getLinearIndex( i, j );
  this->setNode( nodeIdx, x, y );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::getNode(int nodeIdx, double* coordinates) const
{
  assert( coordinates != ATK_NULLPTR );
  assert( nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );

  for ( int i=0; i < this->getDimension(); ++i ) {
    const double* coords = m_coordinates->getCoordinateArray( i );
    coordinates[ i ] = coords[ i ];
  }

}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::getNode( int i, int j, double* coordinates ) const
{
  assert( coordinates != ATK_NULLPTR );
  assert( this->getDimension()==2 );

  const int nodeIdx = m_extent->getLinearIndex( i, j );
  this->getNode( nodeIdx, coordinates );
}

//------------------------------------------------------------------------------
inline
void CurvilinearMesh::getNode(int i, int j, int k, double* coordinates) const
{
  assert( coordinates !=  ATK_NULLPTR );
  assert( this->getDimension()==3 );

  const int nodeIdx = m_extent->getLinearIndex( i, j, k );
  this->getNode( nodeIdx, coordinates );
}

//------------------------------------------------------------------------------
inline double CurvilinearMesh::getNodeCoordinate( int nodeIdx, int idim ) const
{
  assert( nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes() );
  assert( idim >= 0 && idim < this->getDimension() );

  const double* coord = m_coordinates->getCoordinateArray( idim );
  return( coord[nodeIdx] );
}

//------------------------------------------------------------------------------
inline double CurvilinearMesh::getNodeCoordinate( int i, int j, int idim ) const
{
  assert( this->getDimension()==2 );
  assert( idim >= 0 && idim < 2 );

  const int nodeIdx = m_extent->getLinearIndex( i, j );
  return ( this->getNodeCoordinate(nodeIdx, idim) );
}

//------------------------------------------------------------------------------
inline
double CurvilinearMesh::getNodeCoordinate( int i, int j, int k, int idim ) const
{
  assert( this->getDimension()==3 );
  assert( idim >= 0 && idim < 3 );

  const int nodeIdx = m_extent->getLinearIndex( i, j, k );
  return ( this->getNodeCoordinate(nodeIdx, idim) );
}

} /* namespace meshtk */

#endif /* CURVILINEARMESH_HXX_ */
