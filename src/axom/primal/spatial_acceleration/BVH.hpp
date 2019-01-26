/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

#ifndef AXOM_PRIMAL_BVH_H_
#define AXOM_PRIMAL_BVH_H_

#include "axom/config.hpp"        // for Axom compile-time definitions
#include "axom/core/Macros.hpp"   // for Axom macros
#include "axom/core/Types.hpp"    // for fixed bitwidth types

namespace axom
{
namespace primal
{

/*!
 * \brief Enumerates the list of return codes for various BVH operations.
 */
enum BVHReturnCodes
{
  BVH_BUILD_FAILED=-1, //!< indicates that generation of the BVH failed
  BVH_BUILD_OK,        //!< indicates that the BVH was generated successfully
};

/*!
 * \class BVH
 *
 * \brief Defines a Bounding Volume Hierarchy (BVH) spatial acceleration
 *  data structure over a set of geometric entities.
 *
 * The BVH class provides functionality for generating a hierarchical spatial
 * partitioning over a set of geometric entities. Each entity in the BVH is
 * represented by a bounding volume, in this case an axis-aligned bounding box.
 * Once the BVH structure is generated, it is used to accelerate various spatial
 * queries, such as, collision detection, ray tracing, etc., by reducing the
 * search space for a given operation to an abbreviated list of candidate
 * geometric entities to check for a given query.
 *
 */
class BVH
{
public:

  /*!
   * \brief Creates a BVH instance, of specified dimension, over a given set
   *  of geometric entities, each represented by its corresponding axis-aligned
   *  bounding box.
   *
   * \param [in] dimension the spatial dimension for the BVH
   * \param [in] boxes buffer consisting of bounding boxes for each entity.
   * \param [in] numItems the total number of items to store in the BVH.
   *
   * \note boxes is an array of length 2*dimension*numItems, that stores the
   *  two corners of the axis-aligned bounding box corresponding to a given
   *  geometric entity. For example, in 3D, the two corners of the ith bounding
   *  box are given by:
   *  \code
   *    const int offset = i*6;
   *
   *    double xmin = boxes[ offset   ];
   *    double ymin = boxes[ offset+1 ];
   *    double zmin = boxes[ offset+2 ];
   *
   *    double xmax = boxes[ offset+3 ];
   *    double ymax = boxes[ offset+4 ];
   *    double zmax = boxes[ offset+5 ];
   *  \endcode
   *
   * \pre dimension >= 1 && dimension <= 3
   * \pre boxes != nullptr
   * \pre numItems > 0
   */
  BVH( int dimension, const double* boxes, IndexType numItems );

  /*!
   * \brief Destructor.
   */
  ~BVH();

  /*!
   * \brief Generates the BVH
   * \return status set to BVH_BUILD_OK on success.
   */
  int build( );

  /*!
   * \brief Finds the candidate geometric entities that contain a given point.
   *
   * \param [out] candidates buffer consisting of the candidate IDs
   * \param [out] numCandidates the total number of candidates found
   * \param [in] x the x-coordinate of the query point
   * \param [in] y the y-coordinate of the query point
   * \param [in] z the z-coordinate of the query point (needed for 3D only)
   *
   * \note The candidates array is allocated internally by the method and
   *  ownership of the memory is transferred to the caller. Consequently, the
   *  caller is responsible for properly deallocating the candidates buffer.
   *
   * \pre pt != nullptr
   * \pre candidates == nullptr
   */
  void find( IndexType*& candidates,
             IndexType& numCandidates,
             double x,
             double y,
             double z = 0.0f );

  /*!
   * \brief Finds the candidate geometric entities that contain each of the
   *  given query points.
   *
   * \param [out] offsets offset to the candidates array for each query point
   * \param [out] candidates array of the candidate IDs for each query point
   * \param [in]  numPts the total number of query points supplied
   * \param [in]  x array of x-coordinates
   * \param [in]  y array of y-coordinates, may be nullptr if 1D
   * \param [in]  z array of z-coordinates, may be nullptr if 1D or 2D
   *
   * \note the offsets array must be pre-allocated by the caller and
   *
   * \note The candidates array is allocated internally by the method and
   *  ownership of the memory is transferred to the caller. Consequently, the
   *  caller is responsible for properly deallocating the candidates buffer.
   *
   * \pre offsets != nullptr
   * \pre candidates == nullptr
   * \pre x != nullptr
   * \pre y != nullptr if dimension==2 || dimension==3
   * \pre z != nullptr if dimension==3
   */
  void find( IndexType* offsets,
             IndexType*& candidates,
             IndexType numPts,
             const double* x,
             const double* y,
             const double* z=nullptr );

private:

/// \name Private Members
/// @{

  int m_dimension;
  IndexType m_numItems;
  const double* m_boxes;

/// @}

  DISABLE_COPY_AND_ASSIGNMENT(BVH);
  DISABLE_MOVE_AND_ASSIGNMENT(BVH);
};

} /* namespace primal */
} /* namespace axom */

#endif /* AXOM_PRIMAL_BVH_H_ */
