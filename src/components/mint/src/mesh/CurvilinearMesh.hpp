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

  CurvilinearMesh() = delete;

  /*!
   * \brief Constructs a curvilinear mesh instance.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the logical extent of this mesh instance.
   */
  CurvilinearMesh( int dimension, const int64 ext[6] );

  /*!
   * \brief Destructor.
   */
  virtual ~CurvilinearMesh()
  {}

  /// \name SetNode() methods
  /// @{

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \param [in] z the z--coordinate to set at the given node.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes
   * \pre this->getDimension() == 3
   */
  void setNode( IndexType nodeIdx, double x, double y, double z );

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
  void setNode( IndexType i, IndexType j, IndexType k, double x, double y,
                double z );

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \pre this->getDimension() == 2
   */
  void setNode( IndexType nodeIdx, double x, double y );

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] x the x--coordinate to set at the given node.
   * \param [in] y the y--coordinate to set at the given node.
   * \pre this->getDimension() == 2
   */
  void setNode( IndexType i, IndexType j, double x, double y );

  /*!
   * \brief Sets the coordinates of the node at the given index.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] x the x--coordinate to set at the given node.
   * \pre this->getDimension() == 1
   */
  void setNode( IndexType nodeIdx, double x );

  /// @}

  /*!
   * \brief Copy the coordinates of the given node into the provided buffer.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] coords the buffer to copy the coordinates into, of length at
   *  least getDimension().
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != AXOM_NULLPTR
   */
  virtual void getNode( IndexType nodeID, double* coords ) const override final
  { m_coordinates.getCoordinates( nodeID, coords ); }

  /*!
   * \brief Return a pointer to the array of nodal coordinates of the
   *  given dimension.
   *
   * \param [in] dim the dimension to return.
   *
   * \pre 0 <= dim < getDimension()
   */
  virtual double* getCoordinateArray( int dim ) final override
  { return m_coordinates.getCoordinateArray( dim ); }

  /*!
   * \brief Return a constant pointer to the array of nodal coordinates of the
   *  given dimension.
   *
   * \param [in] dim the dimension to return.
   *
   * \pre 0 <= dim < getDimension()
   */
  virtual const double* getCoordinateArray( int dim ) const final override
  { return m_coordinates.getCoordinateArray( dim ); }

private:

  MeshCoordinates m_coordinates;


  CurvilinearMesh(const CurvilinearMesh&); // Not implemented
  CurvilinearMesh& operator=(const CurvilinearMesh&); // Not implemented
};


//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( IndexType nodeIdx, double x, double y,
                                      double z)
{
  SLIC_ASSERT( nodeIdx >= 0 && nodeIdx < getNumberOfNodes() );
  SLIC_ASSERT( getDimension() == 3 );

  m_coordinates.set( nodeIdx, x, y, z );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( IndexType i, IndexType j, IndexType k,
                                      double x, double y, double z )
{
  SLIC_ASSERT( getDimension() == 3 );

  IndexType nodeIdx = m_extent.getLinearIndex( i, j, k );
  setNode( nodeIdx, x, y, z );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( IndexType nodeIdx, double x, double y )
{
  SLIC_ASSERT(  getDimension() == 2 );

  m_coordinates.set( nodeIdx, x, y );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( IndexType i, IndexType j, double x,
                                      double y )
{
  SLIC_ASSERT( getDimension() == 2 );

  IndexType nodeIdx = m_extent.getLinearIndex( i, j );
  setNode( nodeIdx, x, y );
}

//------------------------------------------------------------------------------
inline void CurvilinearMesh::setNode( IndexType nodeIdx, double x )
{
  SLIC_ASSERT( getDimension() == 1 );

  m_coordinates.set( nodeIdx, x );
}

} /* namespace mint */
} /* namespace axom */

#endif /* CURVILINEARMESH_HXX_ */
