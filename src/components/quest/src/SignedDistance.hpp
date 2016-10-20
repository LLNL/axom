/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef SIGNEDDISTANCE_HPP_
#define SIGNEDDISTANCE_HPP_

// ATK includes
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"

#include "quest/BVHTree.hpp"
#include "quest/BoundingBox.hpp"
#include "quest/Orientation.hpp"
#include "quest/Point.hpp"
#include "quest/Triangle.hpp"
#include "quest/Vector.hpp"

#include "mint/Field.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/Mesh.hpp"

#include "common/Utilities.hpp"

// C/C++ includes
#include <cmath> // for std::sqrt()

using namespace asctoolkit;

namespace quest {

template < int NDIMS >
class SignedDistance
{
public:
  typedef Point< double,NDIMS > PointType;
  typedef Vector< double,NDIMS > VectorType;
  typedef Triangle< double,NDIMS > TriangleType;
  typedef BoundingBox< double,NDIMS > BoxType;
  typedef BVHTree< int,NDIMS > BVHTreeType;

public:

  /*!
   *****************************************************************************
   * \brief Creates a SignedDistance instance for queries on the given mesh.
   * \param [in] surfaceMesh user-supplied surface mesh.
   * \param [in] maxObjects max number of objects for spatial decomposition.
   * \param [in] maxLevels max levels for spatial decomposition (optional).
   * \note Default maxLevels is 5 if not specified.
   * \pre surfaceMesh != ATK_NULLPTR
   *****************************************************************************
   */
  SignedDistance( mint::Mesh* surfaceMesh, int maxObjects, int maxLevels=5 );

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  ~SignedDistance();

  /*!
   *****************************************************************************
   * \brief Computes the distance of the given point to the surface mesh.
   * \param [in] queryPnt user-supplied point.
   * \return minDist the signed minimum distance to the surface mesh.
   *****************************************************************************
   */
  double computeDistance( const PointType& queryPnt ) const;

  /*!
   *****************************************************************************
   * \brief Computes the distance of the given point to the surface mesh.
   * \note This is an overloaded method that also returns the BVH buckets and
   *  corresponding triangles used to calculate the signed distance.
   * \param [in]  queryPnt user-supplied point.
   * \param [out] bvh_buckets the buckets of the BVH used to satisfy the query.
   * \param [out] triangles the triangles of the BVH used to satisfy the query.
   * \param [out] my_triangles the triangle used to compute the pseudo-normal.
   * \return minDist the minimum signed distance to the surface mesh.
   *****************************************************************************
   */
  double computeDistance( const PointType& queryPnt,
                          std::vector< int >& bvh_buckets,
                          std::vector< int >& triangles,
                          std::vector< int >& my_triangles,
                          PointType& closest_pt ) const;

  /*!
   *****************************************************************************
   * \brief Returns a const reference to the underlying bucket tree.
   * \return ptr pointer to the underlying bucket tree
   * \post ptr != ATK_NULLPTR
   *****************************************************************************
   */
  const BVHTreeType* getBVHTree( ) const { return m_bvhTree; };

private:

  /*!
   *****************************************************************************
   * \brief Computes the bounding box of the given cell on the surface mesh.
   * \param [in] icell the index of the cell on the surface mesh.
   * \return box bounding box of the cell.
   * \pre m_surfaceMesh != ATK_NULLPTR
   * \pre icell >= 0 && icell < m_surfaceMesh->getMeshNumberOfCells()
   *****************************************************************************
   */
  BoxType getCellBoundingBox( int icell );

  /*!
   *****************************************************************************
   * \brief Searches through the list of candidates surface elements and
   *  computes the minimum squared distance and closest point to a surface
   *  element.
   *
   * \param [in]  pt the query point.
   * \param [in]  candidates list of IDs to the candidate surface elements
   * \param [in]  nelems total number of candidates
   * \param [out] cpt the computed closest point (CPT)
   * \param [out] index corresponding candidate ID of the surface element
   * \param [out] cpt_loc location of the CPT, e.g.node, edge
   * \param [out] surface_elements supplied buffer to store surface elements
   * \param [out] elementIds stores corresponding surface element IDs
   * \param [out] closest_pts supplied buffer to store corresponding CPTs.
   * \param [out] clocs supplied buffer to store corresponding CPT locations.
   *
   * \note The list of candidates is provided by getCandidateSurfaceElements
   *
   * \see SignedDistance::getCandidateSurfaceElements()
   * \see SignedDistance::computeDistance()
   *
   * \pre candidates != ATK_NULLPTR
   * \pre surface_elements != ATK_NULLPTR
   * \pre elementIds != ATK_NULLPTR
   * \pre closest_pts != ATK_NULLPTR
   * \pre clocs != ATK_NULLPTR
   *
   * \post index >= 0 && index < nelems
   *
   * \return dist minimum squared distance to the surface.
   *****************************************************************************
   */
  double getMinSqDistance( const PointType& pt,
                           const int* candidates,
                           int nelems,
                           PointType& cpt,
                           int& index,
                           int& cpt_loc,
                           TriangleType* surface_elements,
                           int* elementIds,
                           PointType* closest_pts,
                           int* clocs ) const;

  /*!
   *****************************************************************************
   * \brief Computes the sign of the given query point given the closest point
   *  and surface element on the surface mesh and neighboring surface elements.
   *
   * \param [in] pt the query point
   * \param [in] cpt the pre-computed closest point (CPT) on the surface
   * \param [in] index the closest surface element candidate list index
   * \param [in] cpt_loc location of CPT on the surface element
   * \param [in] nelems total number of candidates
   * \param [in] elementIds corresponding IDs on the surface mesh
   * \param [in] surface_elements list of surface elements
   * \param [in] closest_pts corresponding list of closest points.
   * \param [in] clocs corresponding location of closest points.
   * \param [out] my_elements list of element used to calculate pseudo-normal
   *
   * \return sign the calculated sign, 1.0 if outside, -1.0 if inside
   *
   * \pre elementIds != ATK_NULLPTR
   * \pre surface_elements != ATK_NULLPTR
   * \pre closest_pts != ATK_NULLPTR
   * \pre clocs != ATK_NULLPTR
   *****************************************************************************
   */
  double computeSign( const PointType& pt,
                      const PointType& cpt,
                      int index,
                      int cpt_loc,
                      int nelems,
                      const int* elementIds,
                      const TriangleType* surface_elements,
                      const PointType* closest_pts,
                      const int* clocs,
                      std::vector< int >& my_elements ) const;

  /*!
   *****************************************************************************
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
   *  \pre bins != ATK_NULLPTR
   *  \pre m_surfaceMesh != ATK_NULLPTR
   *  \pre surface_elements != ATK_NULLPTR
   *  \pre indx != ATK_NULLPTR
   *****************************************************************************
   */
  void getCandidateSurfaceElements( const PointType& pt,
                                    const int* bins,
                                    int nbins,
                                    std::vector< int >& candidates ) const;

  /*!
   *****************************************************************************
   * \brief Returns the maximum distance between a given point and box.
   * \param [in] b user-supplied axis-aligned bounding box.
   * \param [in] pt user-supplied point.
   * \return d maximum distance from a point to a box.
   *****************************************************************************
   */
  double getMaxSqDistance( const BoxType& b, const PointType& pt ) const;

  /*!
   *****************************************************************************
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent its use from the calling application.
   *****************************************************************************
   */
  SignedDistance(): m_surfaceMesh(ATK_NULLPTR), m_bvhTree(ATK_NULLPTR) { };

private:

  mint::Mesh* m_surfaceMesh;     /*!< User-supplied surface mesh. */
  BoxType m_boxDomain;           /*!< bounding box containing surface mesh */
  BVHTreeType* m_bvhTree;        /*!< Spatial acceleration data-structure. */

  DISABLE_COPY_AND_ASSIGNMENT( SignedDistance );

};

} /* namespace quest */

//------------------------------------------------------------------------------
//           SignedDistance Implementation
//------------------------------------------------------------------------------
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

} /* end detail namespace */

//------------------------------------------------------------------------------
template < int NDIMS >
SignedDistance< NDIMS >::SignedDistance(
        mint::Mesh* surfaceMesh, int maxObjects, int maxLevels )
{
  // Sanity checks
  SLIC_ASSERT( surfaceMesh != ATK_NULLPTR );
  SLIC_ASSERT( maxLevels >= 1 );

  m_surfaceMesh    = surfaceMesh;
  const int ncells = m_surfaceMesh->getMeshNumberOfCells();
  const int nnodes = m_surfaceMesh->getMeshNumberOfNodes();

  // compute bounding box of surface mesh
  // NOTE: this should be changed to an oriented bounding box in the future.
  PointType pt;
  for ( int inode=0; inode < nnodes; ++inode ) {
     m_surfaceMesh->getMeshNode( inode, pt.data() );
     m_boxDomain.addPoint( pt );
  }

  // Initialize BucketTree with the surface elements.
  m_bvhTree  = new BVHTreeType( ncells, maxLevels );

  for ( int icell=0; icell < ncells; ++icell ) {
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
  m_bvhTree = ATK_NULLPTR;
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
inline double SignedDistance< NDIMS >::computeDistance( const PointType& pt,
                                          std::vector< int >& buckets,
                                          std::vector< int >& elementIds,
                                          std::vector< int >& my_elements,
                                          PointType& closest_pt ) const
{
  SLIC_ASSERT( m_surfaceMesh != ATK_NULLPTR );
  SLIC_ASSERT( m_bvhTree != ATK_NULLPTR );

  // STEP 0: get list of buckets to satisfy point query
  m_bvhTree->find( pt, buckets );

  // STEP 1: get candidate surface elements
  int nbuckets = buckets.size();
  std::vector< int > candidates;
  this->getCandidateSurfaceElements( pt, &buckets[0], nbuckets, candidates );

  const int nelems = candidates.size();
  elementIds.resize( nelems );

  // STEP 2: process surface elements and compute minimum distance,
  // corresponding closest point and
  int cpt_loc = 0;
  int index   = -1;
  std::vector< TriangleType > triangles;
  std::vector< PointType >    closest_pts;
  std::vector< int >          cpt_locs;
  triangles.resize( nelems );
  closest_pts.resize( nelems );
  cpt_locs.resize( nelems );
  double minSqDist =
      this->getMinSqDistance( pt, &candidates[0], nelems, closest_pt, index,
       cpt_loc, &triangles[0], &elementIds[0], &closest_pts[0], &cpt_locs[0] );

  // STEP 3: compute sign
  double sign =
      this->computeSign( pt, closest_pt, index, cpt_loc, nelems, &elementIds[0],
                &triangles[0], &closest_pts[0], &cpt_locs[0], my_elements );

  // STEP 4: return computed signed distance
  return ( sign*std::sqrt( minSqDist) );
}

//------------------------------------------------------------------------------
template < int NDIMS >
double SignedDistance< NDIMS >::computeSign(
                                       const PointType& pt,
                                       const PointType& cpt,
                                       int index,
                                       int cpt_loc,
                                       int nelems,
                                       const int* elementIds,
                                       const TriangleType* surface_elements,
                                       const PointType* closest_pts,
                                       const int* clocs,
                                       std::vector< int >& my_elements ) const
{
  // Sanity checks
  SLIC_ASSERT( elementIds != ATK_NULLPTR );
  SLIC_ASSERT( surface_elements != ATK_NULLPTR );
  SLIC_ASSERT( closest_pts != ATK_NULLPTR );
  SLIC_ASSERT( clocs != ATK_NULLPTR );

  // STEP 0: if point is outside the bounding box of the surface mesh, then
  // it is outside, just return 1.0
  if ( ! m_boxDomain.contains( pt ) ) {
    /* short-circuit, pt outside bounding box of the surface mesh */
    return 1.0;
  }

  // STEP 1: Otherwise, calculate pseudo-normal N, at the closest point to
  // calculate the sign. There are effectively 3 cases based on the location
  // of the closest point

  VectorType N; // the pseudo-normal, computed below

  if ( cpt_loc == TriangleType::NUM_TRI_VERTS ) {

    // CASE 1: closest point is on the face of the surface element
    N = surface_elements[ index ].normal();
    my_elements.push_back( elementIds[ index ] );

  } else if ( cpt_loc < 0 ) {

    // CASE 2: closest point is on an edge, sum normals of adjacent facets
    for ( int i=0; i < nelems; ++i ) {

      double dist = quest::squared_distance(cpt,closest_pts[i]);
      if ( utilities::isNearlyEqual( dist, 0.0 ) ) {
        N += surface_elements[ i ].normal();
        my_elements.push_back( elementIds[ i ] );

      }
    } // END for

  } else {

    // CASE 3: closest point is on a node, use angle weighted pseudo-normal
    for ( int i=0; i < nelems; ++i ) {

      double dist = quest::squared_distance(cpt,closest_pts[i]);
      if ( utilities::isNearlyEqual( dist, 0.0 ) ) {

        double alpha = surface_elements[ i ].angle( clocs[i] );
        N += ( surface_elements[ i ].normal().unitVector()*alpha );
        my_elements.push_back( elementIds[i] );

      }
    } // END for

  }

  // STEP 2: Given the pseudo-normal, N, and the vector r from the closest point
  // to the query point, compute the sign by checking the sign of their dot
  // product.
  VectorType r( cpt, pt );
  double dotprod = r.dot( N );
  double sign = ( dotprod >= 0.0 )? 1.0 : -1.0;
  SLIC_ASSERT( sign==-1.0 || sign==1.0 );

  return sign;
}

//------------------------------------------------------------------------------
template < int NDIMS >
double SignedDistance< NDIMS >::getMinSqDistance( const PointType& pt,
                                                  const int* candidates,
                                                  int nelems,
                                                  PointType& cpt,
                                                  int& index,
                                                  int& cpt_loc,
                                                  TriangleType* surface_elements,
                                                  int* elementIds,
                                                  PointType* closest_pts,
                                                  int* clocs ) const
{
  SLIC_ASSERT( candidates != ATK_NULLPTR );
  SLIC_ASSERT( surface_elements != ATK_NULLPTR );
  SLIC_ASSERT( closest_pts != ATK_NULLPTR );
  SLIC_ASSERT( clocs != ATK_NULLPTR );

  double minSqDist = std::numeric_limits< double >::max();

  for ( int i=0;  i < nelems; ++i ) {

     const int cellIdx = candidates[ i ];
     elementIds[ i ]   = cellIdx;

     int cellIds[3];
     m_surfaceMesh->getMeshCell( cellIdx, cellIds );

     TriangleType surface_element;
     m_surfaceMesh->getMeshNode( cellIds[0], surface_element[0].data() );
     m_surfaceMesh->getMeshNode( cellIds[1], surface_element[1].data() );
     m_surfaceMesh->getMeshNode( cellIds[2], surface_element[2].data() );
     surface_elements[ i ] = surface_element;

     int iloc = 0;
     const PointType icpt = quest::closest_point( pt, surface_element, &iloc );
     closest_pts[ i ]     = icpt;  // save cpt
     clocs[ i ]           = iloc;  // save cpt location

     double sqDist = quest::squared_distance( pt,icpt );
     if ( sqDist < minSqDist ) {

        minSqDist = sqDist;
        cpt       = icpt;
        cpt_loc   = iloc;
        index     = i;
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
  for ( int i=0; i < npoints; ++i ) {
    dist = std::max( dist, quest::squared_distance(pt,pnts[i]) );
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
  SLIC_ASSERT( m_surfaceMesh != ATK_NULLPTR );
  SLIC_ASSERT( bins != ATK_NULLPTR );

  // STEP 0: count total number of surface elements
  int nelems = 0;
  for ( int ibin=0; ibin < nbins; ++ibin ) {
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
  for ( int ibin=0; ibin < nbins; ++ibin ) {

     const int binIdx   = bins[ ibin ];
     const int nobjects = m_bvhTree->getBucketNumObjects( binIdx );
     const int* objList = m_bvhTree->getBucketObjectArray( binIdx );

     for ( int iobject=0; iobject < nobjects; ++iobject ) {

        const int objectId = objList[ iobject ];
        const int cellIdx  = m_bvhTree->getObjectData( objectId );
        BoxType bbox       = m_bvhTree->getObjectBox( objectId );

//        dist[ icount ]           = this->getMaxSqDistance( bbox, pt );
        dist[ icount ]           = quest::squared_distance( pt, bbox );
        objectIds[ icount ]      = objectId;
        surface_elements[icount] = cellIdx;
        indx[ icount ]           = icount;
        ++icount;

     } // END for all objects within the bin

  } // END for all bins

  SLIC_ASSERT( icount==nelems );

  // STEP 2: sort the indices based on distance
  std::sort( indx, indx+nelems, detail::SortByDistance(dist) );

  // STEP 3: filter out candidates
  candidates.push_back( surface_elements[ indx[0] ] );
  BoxType bbox   = m_bvhTree->getObjectBox( objectIds[ indx[0] ] );
  for ( int i=1;  i < nelems; ++i ) {
    const int idx = indx[i];
    BoxType ib = m_bvhTree->getObjectBox( objectIds[ idx ] );
    if ( bbox.intersects( ib )  ) {
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
inline BoundingBox<double,NDIMS>
SignedDistance< NDIMS >::getCellBoundingBox( int icell )
{
  // Sanity checks
  SLIC_ASSERT( m_surfaceMesh != ATK_NULLPTR );
  SLIC_ASSERT( icell >= 0 && icell < m_surfaceMesh->getMeshNumberOfCells() );

  // Get the cell type, for now we support linear triangle,quad in 3-D and
  // line segments in 2-D.
  const int cellType = m_surfaceMesh->getMeshCellType( icell );
  SLIC_ASSERT( cellType == MINT_TRIANGLE ||
               cellType == MINT_QUAD ||
               cellType == MINT_SEGMENT );
  const int nnodes = mint::cell::num_nodes[ cellType ];

  // Get the cell node IDs that make up the cell
  int* cellIds = new int[ nnodes ];
  m_surfaceMesh->getMeshCell( icell, cellIds );

  // compute the cell's bounding box
  BoxType bb;
  PointType pt;

  for ( int i=0; i < nnodes; ++i ) {

     m_surfaceMesh->getMeshNode( cellIds[ i ], pt.data() );
     bb.addPoint( pt );

  } // END for all cell nodes

  // clean up all dynamically allocated memory
  delete [] cellIds;

  return ( bb );
}

} /* namespace quest */
#endif /* SIGNEDDISTANCE_HPP_ */
