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
   * \brief Constructs a rectilinear mesh instance.
   * \param [in] dimension the dimension of the mesh.
   * \param [in] ext the mesh extent.
   * \param [in] blockId the block ID.
   * \param [in] partitionId the partition ID.
   */
  RectilinearMesh( int dimension, int64 ext[6], int blockId,
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
  inline void setCoordinate( int dim, IndexType i, double coord ) const;

  /*!
   * \brief Gets the coordinate along the given dimension.
   * \return the position of the ith coordinate along the given dimension.
   * \pre dim >= 0 && dim < this->getDimension()
   * \pre i >= 0 && i < ndims[ i ]
   */
  inline double getCoordinate( int dim, IndexType i ) const;

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

} /* namespace mint */
} /* namespace axom */

#endif /* RECTILINEARMESH_HXX_ */
