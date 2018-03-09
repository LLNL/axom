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

#ifndef SIGNEDDISTANCE_HPP_
#define SIGNEDDISTANCE_HPP_

// axom includes
#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "axom_utils/Utilities.hpp"

#include "primal/BVHTree.hpp"
#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/Vector.hpp"

#include "mint/Field.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/Mesh.hpp"

// C/C++ includes
#include <cmath> // for std::sqrt()

namespace axom
{
namespace quest
{

template < int NDIMS >
class SignedDistance
{
public:
  typedef axom::primal::Point< double,NDIMS > PointType;
  typedef axom::primal::Vector< double,NDIMS > VectorType;
  typedef axom::primal::Triangle< double,NDIMS > TriangleType;
  typedef axom::primal::BoundingBox< double,NDIMS > BoxType;
  typedef axom::primal::BVHTree< int,NDIMS > BVHTreeType;

private:

  /// @{
  /// \name Internal Datatype Definitions

  struct cpt_data
  {
    PointType closest_point;
    int candidate_index;
    int cpt_location;

    int nelems;
    std::vector< TriangleType > surface_elements;
    std::vector< int > element_ids;
    std::vector< PointType > closest_pts;
    std::vector< int > cpt_locs;
  };

  /// @}

public:

  /*!
   * \brief Creates a SignedDistance instance for queries on the given mesh.
   * \param [in] surfaceMesh user-supplied surface mesh.
   * \param [in] maxObjects max number of objects for spatial decomposition.
   * \param [in] maxLevels max levels for spatial decomposition (optional).
   * \note Default maxLevels is 5 if not specified.
   * \pre surfaceMesh != AXOM_NULLPTR
   */
  SignedDistance( axom::mint::Mesh* surfaceMesh, int maxObjects,
                  int maxLevels=5 );

  /*!
   * \brief Destructor.
   */
  ~SignedDistance();

  /*!
   * \brief Computes the distance of the given point to the surface mesh.
   * \param [in] queryPnt user-supplied point.
   * \return minDist the signed minimum distance to the surface mesh.
   */
  double computeDistance( const PointType& queryPnt ) const;

  /*!
   * \brief Computes the distance of the given point to the surface mesh.
   * \note This is an overloaded method that also returns the BVH buckets and
   *  corresponding triangles used to calculate the signed distance.
   *
   * \param [in]  queryPnt user-supplied point.
   * \param [out] bvh_buckets the buckets of the BVH used to satisfy the query.
   * \param [out] triangles the triangles of the BVH used to satisfy the query.
   * \param [out] my_triangles the triangle used to compute the pseudo-normal.
   *
   * \note The variables 'triangles'/'my_triangles' are relevant in debug mode.
   *
   * \return minDist the minimum signed distance to the surface mesh.
   */
  double computeDistance( const PointType& queryPnt,
                          std::vector< int >& bvh_buckets,
                          std::vector< int >& triangles,
                          std::vector< int >& my_triangles,
                          PointType& closest_pt ) const;

  /*!
   * \brief Returns a const reference to the underlying bucket tree.
   * \return ptr pointer to the underlying bucket tree
   * \post ptr != AXOM_NULLPTR
   */
  const BVHTreeType* getBVHTree( ) const { return m_bvhTree; };

private:

  /*!
   * \brief Computes the bounding box of the given cell on the surface mesh.
   * \param [in] icell the index of the cell on the surface mesh.
   * \return box bounding box of the cell.
   * \pre m_surfaceMesh != AXOM_NULLPTR
   * \pre icell >= 0 && icell < m_surfaceMesh->getMeshNumberOfCells()
   */
  BoxType getCellBoundingBox( int icell );

  /*!
   * \brief Searches through the list of candidates surface elements and
   *  computes the minimum squared distance and closest point to a surface
   *  element.
   *
   * \param [in]  pt the query point.
   * \param [in]  candidates list of IDs to the candidate surface elements
   * \param [in]  nelems total number of candidates
   * \param [out] cpt_data closest point data
   *
   * \note The list of candidates is provided by getCandidateSurfaceElements
   *
   * \see SignedDistance::getCandidateSurfaceElements()
   * \see SignedDistance::computeDistance()
   *
   * \pre candidates != AXOM_NULLPTR
   * \pre surface_elements != AXOM_NULLPTR
   * \pre elementIds != AXOM_NULLPTR
   * \pre closest_pts != AXOM_NULLPTR
   * \pre clocs != AXOM_NULLPTR
   *
   * \post index >= 0 && index < nelems
   *
   * \return dist minimum squared distance to the surface.
   */
  double getMinSqDistance( const PointType& pt,
                           const int* candidates,
                           int nelems,
                           cpt_data* cpt) const;

  /*!
   * \brief Computes the sign of the given query point given the closest point
   *  and surface element on the surface mesh and neighboring surface elements.
   *
   * \param [in] pt the query point
   * \param [in] cpt pointer to the closest point data for this point query.
   * \param [out] my_elements list of element used to calculate pseudo-normal
   *
   * \return sign the calculated sign, 1.0 if outside, -1.0 if inside
   *
   * \pre elementIds != AXOM_NULLPTR
   * \pre surface_elements != AXOM_NULLPTR
   * \pre closest_pts != AXOM_NULLPTR
   * \pre clocs != AXOM_NULLPTR
   */
  double computeSign( const PointType& pt,
                      const cpt_data* cpt,
                      std::vector< int >& my_elements ) const;

  /*!
   * \brief Returns a sorted list of the candidate surface elements.
   *
   * \param [in] pt the query point.
   * \param [in] bins array of the bins to check.
   * \param [in] nbins size of the bins array.
   * \param [out] surface_elements buffer to store candidate surface elements
   * \param [out] indx indirection array to access surface elements
   * \param [in] nelems total number of elements
   *
   * \note The candidate surface elements are sorted in ascending order, based
   *  on the distance of the query point and the bounding box of the surface
   *  element.
   *
   *  \pre bins != AXOM_NULLPTR
   *  \pre m_surfaceMesh != AXOM_NULLPTR
   *  \pre surface_elements != AXOM_NULLPTR
   *  \pre indx != AXOM_NULLPTR
   */
  void getCandidateSurfaceElements( const PointType& pt,
                                    const int* bins,
                                    int nbins,
                                    std::vector< int >& candidates ) const;

  /*!
   * \brief Returns the maximum distance between a given point and box.
   * \param [in] b user-supplied axis-aligned bounding box.
   * \param [in] pt user-supplied point.
   * \return d maximum distance from a point to a box.
   */
  double getMaxSqDistance( const BoxType& b, const PointType& pt ) const;

  /*!
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent its use from the calling application.
   */
  SignedDistance() : m_surfaceMesh(AXOM_NULLPTR), m_bvhTree(AXOM_NULLPTR) { };

private:
  axom::mint::Mesh* m_surfaceMesh;      /*!< User-supplied surface mesh. */
  BoxType m_boxDomain;           /*!< bounding box containing surface mesh */
  BVHTreeType* m_bvhTree;         /*!< Spatial acceleration data-structure. */

  DISABLE_COPY_AND_ASSIGNMENT( SignedDistance );

};

} // end namespace quest
} // end namespace axom

//------------------------------------------------------------------------------
//           SignedDistance Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace quest
{
namespace detail
{

class SortByDistance
{
public:
  SortByDistance( double* dist ) : m_dist( dist ) { };
  ~SortByDistance() { }
  bool operator()( int i, int j) const { return( m_dist[i] < m_dist[j] ); }
private:
  double* m_dist;
};

} // end namespace detail

//------------------------------------------------------------------------------
template < int NDIMS >
SignedDistance< NDIMS >::SignedDistance(
  axom::mint::Mesh* surfaceMesh, int maxObjects, int maxLevels )
{
  // Sanity checks
  SLIC_ASSERT( surfaceMesh != AXOM_NULLPTR );
  SLIC_ASSERT( maxLevels >= 1 );

  m_surfaceMesh    = surfaceMesh;
  const int ncells = m_surfaceMesh->getMeshNumberOfCells();
  const int nnodes = m_surfaceMesh->getMeshNumberOfNodes();

  // compute bounding box of surface mesh
  // NOTE: this should be changed to an oriented bounding box in the future.
  PointType pt;
  for ( int inode=0 ; inode < nnodes ; ++inode )
  {
    m_surfaceMesh->getMeshNode( inode, pt.data() );
    m_boxDomain.addPoint( pt );
  }

  // Initialize BucketTree with the surface elements.
  m_bvhTree  = new BVHTreeType( ncells, maxLevels );

  for ( int icell=0 ; icell < ncells ; ++icell )
  {
    m_bvhTree->insert( this->getCellBoundingBox( icell ), icell );
  } // END for all cells

  // Build bounding volume hierarchy
  m_bvhTree->build( maxObjects );
}

//------------------------------------------------------------------------------
template < int NDIMS >
SignedDistance< NDIMS >::~SignedDistance( )
{
  delete m_bvhTree;
  m_bvhTree = AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
template < int NDIMS >
inline double SignedDistance< NDIMS >::computeDistance(
  const PointType& pt) const
{
  std::vector< int > buckets;
  std::vector< int > elements;
  std::vector< int > my_elements;
  PointType closest_pt;
  double dist = this->computeDistance( pt, buckets, elements, my_elements,
                                       closest_pt );
  return ( dist );
}

//------------------------------------------------------------------------------
template < int NDIMS >
inline double SignedDistance< NDIMS >::computeDistance(
  const PointType& pt,
  std::vector< int >& buckets,
  std::vector< int >& AXOM_DEBUG_PARAM(elementIds),
  std::vector< int >& my_elements,
  PointType& closest_pt ) const
{
  SLIC_ASSERT( m_surfaceMesh != AXOM_NULLPTR );
  SLIC_ASSERT( m_bvhTree != AXOM_NULLPTR );

  // STEP 0: get list of buckets to satisfy point query
  m_bvhTree->find( pt, buckets );

  // STEP 1: get candidate surface elements
  int nbuckets = buckets.size();
  std::vector< int > candidates;
  this->getCandidateSurfaceElements( pt, &buckets[0], nbuckets, candidates );

  const int nelems = candidates.size();

  // STEP 2: process surface elements and compute minimum distance and
  // corresponding closest point.
  cpt_data cpt;
  double minSqDist = this->getMinSqDistance( pt, &candidates[0], nelems, &cpt );
  closest_pt = cpt.closest_point;
#ifdef AXOM_DEBUG
  elementIds = cpt.element_ids;
#endif

  // STEP 3: compute sign
  double sign = this->computeSign( pt, &cpt, my_elements );

  // STEP 4: return computed signed distance
  return ( sign*std::sqrt( minSqDist) );
}

//------------------------------------------------------------------------------
template < int NDIMS >
double SignedDistance< NDIMS >::computeSign(
  const PointType& pt,
  const cpt_data* cpt,
  std::vector< int >& AXOM_DEBUG_PARAM(my_elements) ) const
{
  // Sanity checks
  SLIC_ASSERT( cpt != AXOM_NULLPTR );

  // STEP 0: if point is outside the bounding box of the surface mesh, then
  // it is outside, just return 1.0
  if ( !m_boxDomain.contains( pt ) )
  {
    /* short-circuit, pt outside bounding box of the surface mesh */
    return 1.0;
  }

  // STEP 1: Otherwise, calculate pseudo-normal N, at the closest point to
  // calculate the sign. There are effectively 3 cases based on the location
  // of the closest point

  VectorType N; // the pseudo-normal, computed below

  const int cpt_loc = cpt->cpt_location;
  const int index   = cpt->candidate_index;
  const int nelems  = cpt->nelems;

  if ( cpt_loc >= TriangleType::NUM_TRI_VERTS )
  {

    // CASE 1: closest point is on the face of the surface element
    N = cpt->surface_elements[ index ].normal();
#ifdef AXOM_DEBUG
    my_elements.push_back( cpt->element_ids[ index ] );
#endif

  }
  else if ( cpt_loc < 0 )
  {

    // CASE 2: closest point is on an edge, sum normals of adjacent facets
    for ( int i=0 ; i < nelems ; ++i )
    {

      double dist = axom::primal::squared_distance( cpt->closest_point,
                                                    cpt->closest_pts[i]  );

      if ( utilities::isNearlyEqual( dist, 0.0 ) )
      {
        N += cpt->surface_elements[ i ].normal();
#ifdef AXOM_DEBUG
        my_elements.push_back( cpt->element_ids[ i ] );
#endif
      } // END if

    } // END for

  }
  else
  {

    // CASE 3: closest point is on a node, use angle weighted pseudo-normal
    for ( int i=0 ; i < nelems ; ++i )
    {

      double dist = axom::primal::squared_distance( cpt->closest_point,
                                                    cpt->closest_pts[i] );

      if ( utilities::isNearlyEqual( dist, 0.0 ) )
      {
        double alpha = cpt->surface_elements[ i ].angle( cpt->cpt_locs[ i ] );
        N += ( cpt->surface_elements[ i ].normal().unitVector()*alpha );
#ifdef AXOM_DEBUG
        my_elements.push_back( cpt->element_ids[ i ] );
#endif
      } // END if

    } // END for

  }

  // STEP 2: Given the pseudo-normal, N, and the vector r from the closest point
  // to the query point, compute the sign by checking the sign of their dot
  // product.
  VectorType r( cpt->closest_point, pt );
  double dotprod = r.dot( N );
  double sign = ( dotprod >= 0.0 ) ? 1.0 : -1.0;
  SLIC_ASSERT( sign==-1.0 || sign==1.0 );

  return sign;
}

//------------------------------------------------------------------------------
template < int NDIMS >
double SignedDistance< NDIMS >::getMinSqDistance( const PointType& pt,
                                                  const int* candidates,
                                                  int nelems,
                                                  cpt_data* cpt ) const
{
  SLIC_ASSERT( candidates != AXOM_NULLPTR );
  SLIC_ASSERT( cpt != AXOM_NULLPTR );

  cpt->nelems = nelems;
  cpt->surface_elements.resize( nelems );
  cpt->closest_pts.resize( nelems );
  cpt->cpt_locs.resize( nelems );
  cpt->element_ids.resize( nelems );

  PointType* closest_pts = &(cpt->closest_pts)[0];
  int* cpt_locs          = &(cpt->cpt_locs)[0];

  double minSqDist = std::numeric_limits< double >::max();

  for ( int i=0 ; i < nelems ; ++i )
  {

    const int cellIdx = candidates[ i ];
    cpt->element_ids[ i ] = cellIdx;

    int cellIds[3];
    m_surfaceMesh->getMeshCell( cellIdx, cellIds );

    TriangleType surface_element;
    m_surfaceMesh->getMeshNode( cellIds[0], surface_element[0].data() );
    m_surfaceMesh->getMeshNode( cellIds[1], surface_element[1].data() );
    m_surfaceMesh->getMeshNode( cellIds[2], surface_element[2].data() );
    cpt->surface_elements[ i ] = surface_element;

    closest_pts[ i ] =
      axom::primal::closest_point( pt,surface_element,&cpt_locs[i] );

    double sqDist = axom::primal::squared_distance( pt,closest_pts[i] );
    if ( sqDist < minSqDist )
    {

      minSqDist             = sqDist;
      cpt->closest_point    = closest_pts[i];
      cpt->cpt_location     = cpt_locs[i];
      cpt->candidate_index  = i;
    }

  } // END for all elements

  return minSqDist;
}

//------------------------------------------------------------------------------
template < int NDIMS >
double SignedDistance< NDIMS >::getMaxSqDistance( const BoxType& b,
                                                  const PointType& pt ) const
{
  std::vector< PointType > pnts;
  BoxType::getPoints( b, pnts );

  const int npoints = pnts.size();
  SLIC_ASSERT( npoints==4 || npoints==8 );

  double dist = std::numeric_limits< double >::min();
  for ( int i=0 ; i < npoints ; ++i )
  {
    dist = std::max( dist, axom::primal::squared_distance(pt,pnts[i]) );
  } // END for all points

  return dist;
}

//------------------------------------------------------------------------------
template < int NDIMS >
void SignedDistance< NDIMS >::getCandidateSurfaceElements(
  const PointType& pt,
  const int* bins,
  int nbins,
  std::vector< int >& candidates ) const
{
  // Sanity checks
  SLIC_ASSERT( m_surfaceMesh != AXOM_NULLPTR );
  SLIC_ASSERT( bins != AXOM_NULLPTR );

  // STEP 0: count total number of surface elements
  int nelems = 0;
  for ( int ibin=0 ; ibin < nbins ; ++ibin )
  {
    const int bucketIdx = bins[ ibin ];
    nelems += m_bvhTree->getBucketNumObjects( bucketIdx );
  }

  // STEP 0: allocate array to store the distance of the bounding box of each
  // surface element to the query point.
  int* surface_elements = new int[ nelems ];
  double* dist          = new double[ nelems ];
  int* objectIds        = new int[ nelems ];
  int* indx             = new int[ nelems ];

  // STEP 1: get flat array of the surface elements and compute distances
  int icount = 0;
  for ( int ibin=0 ; ibin < nbins ; ++ibin )
  {

    const int binIdx   = bins[ ibin ];
    const int nobjects = m_bvhTree->getBucketNumObjects( binIdx );
    const int* objList = m_bvhTree->getBucketObjectArray( binIdx );

    for ( int iobject=0 ; iobject < nobjects ; ++iobject )
    {

      const int objectId = objList[ iobject ];
      const int cellIdx  = m_bvhTree->getObjectData( objectId );
      BoxType bbox       = m_bvhTree->getObjectBox( objectId );

//        dist[ icount ]           = this->getMaxSqDistance( bbox, pt );
      dist[ icount ]           = axom::primal::squared_distance( pt, bbox );
      objectIds[ icount ]      = objectId;
      surface_elements[icount] = cellIdx;
      indx[ icount ]           = icount;
      ++icount;

    }  // END for all objects within the bin

  } // END for all bins

  SLIC_ASSERT( icount==nelems );

  // STEP 2: sort the indices based on distance
  std::sort( indx, indx+nelems, detail::SortByDistance(dist) );

  // STEP 3: filter out candidates
  candidates.push_back( surface_elements[ indx[0] ] );
  const BoxType& bbox = m_bvhTree->getObjectBox( objectIds[ indx[0] ] );
  double maxDist = this->getMaxSqDistance( bbox, pt );

  for ( int i=1 ; i < nelems ; ++i )
  {
    const int idx = indx[i];
    if ( dist[ idx ] < maxDist )
    {
      candidates.push_back( surface_elements[ idx ] );
    }
  }

  // STEP 4: delete temporary distance array
  delete [] dist;
  delete [] objectIds;
  delete [] indx;
  delete [] surface_elements;
}

//------------------------------------------------------------------------------
template < int NDIMS >
inline axom::primal::BoundingBox< double,NDIMS >
SignedDistance< NDIMS >::getCellBoundingBox( int icell )
{
  // Sanity checks
  SLIC_ASSERT( m_surfaceMesh != AXOM_NULLPTR );
  SLIC_ASSERT( icell >= 0 && icell < m_surfaceMesh->getMeshNumberOfCells() );

  // Get the cell type, for now we support linear triangle,quad in 3-D and
  // line segments in 2-D.
  const int cellType = m_surfaceMesh->getMeshCellType( icell );
  SLIC_ASSERT( cellType == MINT_TRIANGLE ||
               cellType == MINT_QUAD ||
               cellType == MINT_SEGMENT );
  const int nnodes = axom::mint::cell::num_nodes[ cellType ];

  // Get the cell node IDs that make up the cell
  int* cellIds = new int[ nnodes ];
  m_surfaceMesh->getMeshCell( icell, cellIds );

  // compute the cell's bounding box
  BoxType bb;
  PointType pt;

  for ( int i=0 ; i < nnodes ; ++i )
  {

    m_surfaceMesh->getMeshNode( cellIds[ i ], pt.data() );
    bb.addPoint( pt );

  } // END for all cell nodes

  // clean up all dynamically allocated memory
  delete [] cellIds;

  return ( bb );
}

} // end namespace quest
} // end namespace axom

#endif /* SIGNEDDISTANCE_HPP_ */
