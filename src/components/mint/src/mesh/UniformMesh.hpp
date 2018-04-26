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

#ifndef UNIFORMMESH_HXX_
#define UNIFORMMESH_HXX_

#include "mint/StructuredMesh.hpp"        /* for StructuredMesh */
#include "mint/config.hpp"             /* for IndexType, int64 */
#include "slic/slic.hpp"                  /* for SLIC macros */

#include <algorithm>                      /* for std::fill */

namespace axom
{
namespace mint
{

class UniformMesh : public StructuredMesh
{
public:

  /*!
   * \brief Default constructor. Disabled.
   */
  UniformMesh();

  /*!
   * \brief Constructs a uniform mesh defined by the origin, spacing and extent.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] origin the origin coordinates of the mesh
   * \param [in] h the spacing in each dimension.
   * \param [in] ext the extent of this mesh instance.
   */
  UniformMesh( int dimension, const double origin[3], const double h[3],
               const int64 ext[6] );

  /*!
   * \brief Constructs a uniform mesh defined by the extent of the mesh in each
   *  direction and an axis aligned bounding box.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the extent of this mesh instance.
   * \param [in] lower_bound a corner of the bounding box.
   * \param [in] upper_bound the corner opposite lower_bound of the bounding
   * box.
   */
  UniformMesh( int dimension, const int64 ext[6],
               const double lower_bound[3], const double upper_bound[3] );

/// \name Virtual methods
/// @{

  /*!
   * \brief Destructor.
   */
  virtual ~UniformMesh()
  {}

/// \name Nodes
/// @{

  /*!
   * \brief Copy the coordinates of the given node into the provided buffer.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] coords the buffer to copy the coordinates into, of length at
   *  least getDimension().
   *
   * \note provided only for convenience, do not use inside a loop. Instead use
   *  getOrigin() and getSpacing() to calculate the nodal coordinates.
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != AXOM_NULLPTR
   */
  virtual void getNode( IndexType nodeID, double* node ) const override final;

  /*!
   * \brief Return a pointer to the nodal positions in the specified dimension.
   *  Since the UniformMesh holds no such array it return AXOM_NULLPTR.
   */
  /// @{

  virtual double* getCoordinateArray( int AXOM_NOT_USED(dim) ) final override
  { return AXOM_NULLPTR; }

  virtual const double* getCoordinateArray( int AXOM_NOT_USED(dim) )
                                                            const final override
  { return AXOM_NULLPTR; }

  /// @}

/// @}

/// @}

/// \name Attribute Querying Methods
/// @{

  /*!
   * \brief Returns the origin of the Uniform Mesh
   * \param [out] origin user-supplied buffer to store the mesh origin.
   */
  void getOrigin( double origin[3] ) const;

  /*!
   * \brief Returns the spacing of the Uniform Mesh.
   * \param [out] h user-supplied buffer to store the spacing of the mesh.
   */
  void getSpacing( double h[3] ) const;

/// @}

private:

  double m_origin[3];
  double m_h[3];

  DISABLE_COPY_AND_ASSIGNMENT(UniformMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(UniformMesh);
};


//------------------------------------------------------------------------------
//          In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void UniformMesh::getOrigin( double origin[3] ) const
{
  SLIC_ASSERT( origin != AXOM_NULLPTR );
  memcpy( origin, m_origin, 3*sizeof(double) );
}

//------------------------------------------------------------------------------
inline void UniformMesh::getSpacing( double h[3] ) const
{
  SLIC_ASSERT( h != AXOM_NULLPTR );
  memcpy( h, m_h, 3*sizeof(double) );
}

//------------------------------------------------------------------------------
inline void UniformMesh::getNode( IndexType nodeID, double* node ) const
{
  SLIC_ASSERT( 0 <= nodeID && nodeID < getNumberOfNodes() );
  SLIC_ASSERT( node != AXOM_NULLPTR );
  if ( getDimension() == 1 )
  {
    node[0] = m_origin[0] + nodeID * m_h[0];
  }
  else if ( getDimension() == 2 )
  {
    IndexType i, j;
    m_extent->getGridIndex( nodeID, i, j );
    node[0] = m_origin[0] + i * m_h[0];
    node[1] = m_origin[1] + j * m_h[1];
  }
  else
  {
    SLIC_ASSERT( getDimension() == 3 );
    IndexType i, j, k;
    m_extent->getGridIndex( nodeID, i, j, k );
    node[0] = m_origin[0] + i * m_h[0];
    node[1] = m_origin[1] + j * m_h[1];
    node[2] = m_origin[2] + k * m_h[2];
  }
}

} /* namespace mint */
} /* namespace axom */

#endif /* UNIFORMMESH_HXX_ */
