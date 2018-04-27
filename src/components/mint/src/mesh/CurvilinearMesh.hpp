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
   * \brief Default constructor. Disabled.
   */
  CurvilinearMesh() = delete;

  /*!
   * \brief Constructs a curvilinear mesh instance.
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the logical extent of this mesh instance.
   */
  CurvilinearMesh( int dimension, const int64 ext[6] );

/// \name Virtual methods
/// @{

  /*!
   * \brief Destructor.
   */
  virtual ~CurvilinearMesh()
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
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != AXOM_NULLPTR
   */
  virtual void getNode( IndexType nodeID, double* coords ) const override final
  { m_coordinates.getCoordinates( nodeID, coords ); }

  /*!
   * \brief Returns pointer to the requested mesh coordinate buffer.
   *
   * \param [in] dim the dimension of the requested coordinate buffer
   * \return ptr pointer to the coordinate buffer.
   *
   * \note The length of the returned buffer is getNumberOfNodes().
   *
   * \pre dim >= 0 && dim < dimension()
   * \pre dim == X_COORDINATE || dim == Y_COORDINATE || dim == Z_COORDINATE
   */
  /// @{

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

  /// @}

/// @}

/// @}

private:

  MeshCoordinates m_coordinates;

  CurvilinearMesh(const CurvilinearMesh&); // Not implemented
  CurvilinearMesh& operator=(const CurvilinearMesh&); // Not implemented
};

} /* namespace mint */
} /* namespace axom */

#endif /* CURVILINEARMESH_HXX_ */
