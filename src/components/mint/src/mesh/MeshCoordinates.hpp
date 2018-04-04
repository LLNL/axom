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

#ifndef MINT_MESHCOORDINATES_HPP_
#define MINT_MESHCOORDINATES_HPP_

#include "axom/Macros.hpp" // for Axom macros and definitions

#include "mint/config.hpp" // for mint::IndexType
#include "mint/Array.hpp"  // for mint::Array

namespace axom
{

#ifdef MINT_USE_SIDRE

// Sidre forward declarations
namespace sidre
{
class Group;
}

#endif

namespace mint
{

static constexpr int X_COORDINATE = 0;
static constexpr int Y_COORDINATE = 1;
static constexpr int Z_COORDINATE = 2;

/*!
 * \class MeshCoordinates
 *
 * \brief Provides functionality to store and operate on mesh node coordinates.
 *
 *  A MeshCoordinates object is generally tied to an associated Mesh instance
 *  and provides the means to store, access and modify the node coordinates of
 *  a mesh.
 *
 *  The MeshCoordinates object may be constructed using (a) native storage,
 *  (b) external storage, or, (c) Sidre:
 *
 *  * <b> Native Storage </b> <br />
 *
 *    When using native storage, the MeshCoordinates object owns all associated
 *    memory. The storage can dynamically grow as needed, e.g., when adding
 *    more nodes. Typically, extra space is allocated to minimize the number
 *    of re-allocations. At any given instance, the total node capacity
 *    can be queried by calling the capacity() function. The extra memory can
 *    be returned to the system by calling the shrink() method.
 *
 *    When all extra memory is exhausted, appending a new node triggers a
 *    re-allocation. The amount of extra space that is allocated is according
 *    to the <em> resize_ratio </em> parameter, which is set to 2.0 by default.
 *    The <em> resize_ratio </em> may be queried and set to a different value
 *    by the getResizeRatio() and setResizeRatio() functions respectively.
 *
 *    When the MeshCoordinates object goes out-of-scope, all memory associated
 *    with the given instance is returned to the system.
 *
 *  * <b> External Storage </b> <br />
 *
 *    A MeshCoordinates object may also be constructed from external,
 *    user-supplied buffers that store the mesh node coordinates. In this case,
 *    the memory is owned by the caller. The MeshCoordinates object just keeps
 *    pointers to the user-supplied buffers.
 *
 *    \warning Since the memory is not owned by the MeshCoordinates object
 *     when external buffers are supplied, the MeshCoordinates object cannot
 *     dynamically grow the storage. Consequently, the number of nodes the
 *     MeshCoordinates instance can hold is fixed. All calls to `shrink()` and 
 *     `reserve()` will fail.
 *
 *    \warning Moreover, when the MeshCoordinates object goes out-of-scope, the
 *     associated buffers are not deleted. The caller owns the external data
 *     and has the responsibility of properly de-allocating the associated
 *     memory.
 *
 *  * <b> Sidre </b> <br />
 *
 *    A MeshCoordinates object may also be constructed from a sidre::Group which
 *    conforms to the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
 *    mesh blueprint </a>.
 *
 *    A MeshCoordinates object that is bound to a particular sidre::Group
 *    supports all operations transparently including dynamically growing the
 *    storage to hold more nodes as needed, but, instead, Sidre owns the memory.
 *    All memory management operations are delegated to Sidre.
 *
 *    \warning Once the MeshCoordinates object goes out-of-scope, the data
 *     stays remains persistent in Sidre.
 *
 * \warning Reallocations tend to be costly operations in terms of performance.
 *  Use `reserve()` when the number of nodes is known a priori, or opt to
 *  use a constructor that takes an actual size and capacity when possible.
 *
 * \see mint::Array
 * \see mint::Mesh
 * \see sidre::Group
 */
class MeshCoordinates
{
public:

/// \name Native Storage Constructors
/// @{

  /*!
   * \brief Creates an empty MeshCoordinates instance of specified dimension.
   *
   * \param [in] dimension the mesh dimension, i.e., 1, 2 or 3
   *
   * \note The resulting MeshCoordinates object is empty, but, it has initial
   *  default capacity based on an internal settings.
   *
   * \pre 1 <= dimension <= 3
   *
   * \post numNodes() == 0
   * \post empty() == true
   * \post capacity() > 0
   * \post getCoordinateArray( i ) != AXOM_NULLPTR \f$ \forall i \in [0,N-1]\f$
   *  where N is the dimension of the MeshCoordinates object.
   */
  explicit MeshCoordinates( int dimension );

  /*!
   * \brief Creates a MeshCoordinates instance of specified dimension and
   *  given number of nodes. Optionally, an initial max capacity may be
   *  specified in the constructor to reserve additional space for storing
   *  nodes. If a capacity is not specified, an internal initial default
   *  capacity will be computed instead.
   *
   * \param [in] dimension the mesh dimension.
   * \param [in] numNodes the number of nodes this instance will hold initially
   * \param [in] capacity initial max capacity to reserve space for (optional).
   *
   * \pre 1 <= dimension() <= 3
   * \post numNodes() <= capacity()
   * \post numNodes() == numNodes
   * \post if capacity == USE_DEFAULT then
   *  capacity() == max(DEFAULT_CAPACITY, numNodes()*DEFAULT_RESIZE_RATIO)
   */
  MeshCoordinates( int dimension, IndexType numNodes,
                   IndexType capacity=USE_DEFAULT );

/// @}

/// \name External Storage Constructors
/// @{

  /*!
   * \brief Creates a MeshCoordinates object from the given coordinate buffers.
   *
   * \param [in] x pointer to the x-coordinates
   * \param [in] y pointer to the y-coordinates, may be AXOM_NULLPTR if 1D
   * \param [in] z pointer to the z-coordinates, may be AXOM_NULLPTR if 2D
   * \param [in] numNodes the number of nodes in the supplied array
   *
   * \warning This constructor wraps the given coordinate buffers and does not
   *  own the data. Consequently, new nodes cannot be added to this
   *  MeshCoordinates instance when using this constructor.
   *
   * \warning All calls to shrink(), append(), resize() and reserve() will fail.
   *
   * \pre x != AXOM_NULLPTR
   * \pre y != AXOM_NULLPTR if dimension==2 || dimension==3
   * \pre z != AXOM_NULLPTR if dimension==3
   * \pre numNodes >= 1
   *
   * \post 1 <= dimension() <= 3
   * \post empty() == false
   * \post numNodes() == numNodes
   * \post numNodes() <= capacity()
   */
  MeshCoordinates( IndexType numNodes, double* x, double* y=AXOM_NULLPTR,
                   double* z=AXOM_NULLPTR );
/// @}

/// \name Sidre Constructors
/// @{

#ifdef MINT_USE_SIDRE

  /*!
   * \brief Creates a MeshCoordinates object from the given sidre::Group.
   *
   * \param [in] group the sidre::Group instance that holds mesh coordinates.
   *
   * \note The supplied sidre::Group must conform to the mesh blueprint
   *  computational mesh protocol for storing coordinates. For details on
   *  the mesh blueprint, see http://llnl-conduit.readthedocs.io/en/latest/.
   *
   * \pre group != AXOM_NULLPTR
   *
   * \see sidre::Group
   */
  MeshCoordinates( sidre::Group* group );

  /*!
   * \brief Creates a MeshCoordinates object on the given sidre::Group
   *
   * \param [in,out] group the sidre::Group to hold the mesh coordinates
   * \param [in] dimension the mesh dimension
   * \param [in] numNodes the number of nodes initially this instance will hold.
   * \param [in] capacity the initial max node capacity (optional)
   *
   * \note This constructor will create the appropriate groups and views under
   *  the given sidre::Group.
   *
   * \warning The user-supplied sidre::Group is expected to be empty
   *
   * \pre group != AXOM_NULLPTR.
   * \pre group->getNumGroups()==0
   * \pre group->getNumViews()==0
   *
   * \post group->getNumGroups()==2
   * \post group->getNumViews()==dimension
   * \post numNodes() == numNodes
   * \post numNodes() <= capacity()
   *
   * \see sidre::Group
   */
  MeshCoordinates( sidre::Group* group,
                   int dimension,
                   IndexType numNodes,
                   IndexType capacity=USE_DEFAULT );

#endif

/// @}

  /*!
   * \brief Destructor, free's the allocated vectors.
   */
  ~MeshCoordinates();

/// \name Attribute Query Methods
/// @{

  /*!
   * \brief Returns the dimension of this MeshCoordinates instance.
   * \return dim the dimension
   * \post 1 <= dim <= 3
   */
  inline int dimension() const { return m_ndims; };

  /*!
   * \brief Get the maximum number of points that can currently be held.
   * \return N the capacity of m_coordinates.
   * \post N >= numNodes()
   */
  inline IndexType capacity() const { return m_coordinates[0]->capacity(); }

  /*!
   * \brief Returns the number of nodes in this MeshCoordinates instance.
   * \return N the number of nodes in this MeshCoordinates instance.
   */
  inline IndexType numNodes() const { return m_coordinates[0]->size(); }

  /*
   * \brief Checks if this MeshCoordinates instance is empty.
   * \return status true if numNodes()==0, else, false.
   */
  inline bool empty() const { return (numNodes() == 0); }

/// @}

/// \name Resize Ratio Get/Set Methods
/// @{

  /*!
   * \brief Returns the resize ratio by which the capacity will increase.
   * \return N the ratio by which the capacity will increase
   * \see mint::Array::getResizeRatio()
   */
  inline double getResizeRatio() const
  { return m_coordinates[0]->getResizeRatio(); }

  /*!
   * \brief Sets the resize ratio by which the capacity will increase when
   *  a dynamically re-allocation is triggered.
   *
   * \param [in] ratio the ratio by which the capacity will increase
   * \see mint::Array::setResizeRatio()
   */
  inline void setResizeRatio( double ratio );

/// @}

/// \name Modify Methods
/// @{

  /*!
   * \brief Appends a new node to the MeshCoordinates instance
   *
   * \param [in] x pointer to buffer consisting the coordinates of the new node.
   * \return idx the index of the new node
   *
   * \note x should point to a buffer that is at least m_ndims long
   *
   * \warning The append() operation is invalid when a MeshCoordinates object
   *  is constructed using external buffers. Since, in this case the object does
   *  not own the associated memory, dynamic resizing is not allowed.
   *
   * \pre x != AXOM_NULLPTR
   * \post idx >=0 && idx == numNodes()-1
   * \post the number of nodes is incremented by one.
   */
  inline IndexType append( const double* x );

  /*!
   * \brief Sets the coordinates for the given node.
   *
   * \param [in] nodeIdx the index of the node whose coordinates to set
   * \param [in] x pointer to the coordinates buffer.
   *
   * \note x should point to a buffer that is at least m_ndims long
   */
  inline void set( IndexType nodeIdx, const double* x );

/// @}

/// \name Access Methods
/// @{

  /*!
   * \brief Returns the coordinate of a point at the given dimension.
   *
   * \param [in] nodeIdx the index of the point in query.
   * \param [in] dim the dimension in query.
   *
   * \return coord the coordinate of the point.
   *
   * \pre dim < m_ndims
   * \pre (nodeIdx >= 0) && (nodeIdx < this->getNumberOfPoints())
   */
  inline double getCoordinate( IndexType nodeIdx, int dim );

  /*!
   * \brief Returns pointer to the coordinate array at the requested dimension.
   * \param [in] dim the requested dimension.
   * \return coord_array pointer to the coordinate array
   *
   * \pre dim < m_ndims
   * \post coord_array != AXOM_NULLPTR.
   */
  /// @{

  inline double* getCoordinateArray( int dim )
  {
    SLIC_ASSERT( indexInRange( dim, 0, m_ndims-1) );
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
    return m_coordinates[ dim ]->getData();
  }

  inline const double* getCoordinateArray( int dim ) const
  {
    SLIC_ASSERT( indexInRange( dim, 0, m_ndims-1)  );
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
    return m_coordinates[ dim ]->getData();
  }

  /// @}

/// @}

/// \name Size/Capacity Modifiers
/// @{

  /*!
   * \brief Changes the max capacity to the specified capacity.
   *
   * \warning The reserve() operation is invalid when a MeshCoordinates object
   *  is constructed using external buffers. Since, in this case the object does
   *  not own the associated memory, dynamic resizing is not allowed.
   *
   * \param [in] capacity the new max node capacity.
   *
   * \post capacity() == capacity
   * \post numNodes() <= capacity()
   */
  inline void reserve( IndexType capacity );

  /*!
   * \brief Changes the number of nodes to the specified size.
   *
   * \warning The resize() operation is invalid when a MeshCoordinates object
   *  is constructed using external buffers. Since, in this case the object does
   *  not own the associated memory, dynamic resizing is not allowed.
   *
   * \post numNodes()==size
   * \post numNodes() <= capacity()
   */
  inline void resize( IndexType size );

  /*!
   * \brief Returns all extra memory to the system.
   *
   * \warning The shrink() operation is invalid when a MeshCoordinates object
   *  is constructed using external buffers. Since, in this case the object does
   *  not own the associated memory, dynamic resizing is not allowed.
   *
   * \post numNodes() == capacity()
   */
  inline void shrink( );

/// @}

private:

/// \name Private helper methods
/// @{

  /*!
   * \brief Checks if the given index is within the specified range.
   *
   * \param [in] i the index to check
   * \param [in] A lower bound
   * \param [in] B upper bound
   *
   * \return status true \f$ \iff i \in [A,B] \f$, otherwise false
   */
  inline bool indexInRange( int i, int A, int B ) const
  { return ( i >= A && i <= B ); }

  /*!
   * \brief Helper method to check to validate the supplied dimension.
   * \note Used primarily to do pre-condition checks
   * \return status true if the dimension is invalid, otherwise, false.
   */
  inline bool invalidDimension( ) const
  { return (m_ndims < 0 || m_ndims > 3); }

  /*!
   * \brief Checks if the given node index is within the valid range.
   * \param [in] idx the node index to check
   * \return status true if the index is valid, false, otherwise.
   */
  inline bool validIndex( IndexType idx ) const
  { return indexInRange( idx, 0, numNodes() - 1 ); }

  /*!
   * \brief Helper method to initialize the internal array data-structures.
   *
   * \param [in] numNodes the number of nodes to store
   * \param [in] maxCapacity initial max capacity of the arrays
   */
  inline void initialize( IndexType numNodes, IndexType maxCapacity );

  /*!
   * \brief Helper method to check if the internal arrays are consistent.
   * \return status true if the arrays are consistent, else, false.
   */
  bool consistencyCheck() const;

/// @}

  int m_ndims;
  Array< double >* m_coordinates[3] = {AXOM_NULLPTR, AXOM_NULLPTR, AXOM_NULLPTR};

  DISABLE_COPY_AND_ASSIGNMENT( MeshCoordinates );
  DISABLE_MOVE_AND_ASSIGNMENT( MeshCoordinates );
};

//------------------------------------------------------------------------------
// In-lined method implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline IndexType MeshCoordinates::append( const double* x )
{
  SLIC_ASSERT( x != AXOM_NULLPTR );

  IndexType idx = numNodes();

  for ( int dim = 0; dim < m_ndims; ++dim )
  {
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
    m_coordinates[ dim ]->append( x[ dim ] );
  }

  SLIC_ASSERT( idx == numNodes() - 1 );
  SLIC_ASSERT( this->validIndex( idx ) ) ;
  SLIC_ASSERT( this->consistencyCheck( ) );

  return idx;
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::set( IndexType nodeIdx, const double* x )
{
  SLIC_ASSERT( x != AXOM_NULLPTR );
  SLIC_ASSERT( this->validIndex( nodeIdx ) );

  for ( int dim=0; dim < m_ndims; ++dim )
  {
    double* data = m_coordinates[ dim ]->getData( );
    SLIC_ASSERT( m_coordinates[ dim ]->numComponents() == 1 );
    SLIC_ASSERT( data != AXOM_NULLPTR );

    data[ nodeIdx ] = x[ dim ];
  }
}

//------------------------------------------------------------------------------
inline double MeshCoordinates::getCoordinate( IndexType nodeIdx, int dim )
{
  SLIC_ASSERT( dim < m_ndims );
  SLIC_ASSERT( ( nodeIdx >= 0 ) && ( nodeIdx < numNodes() ) );
  SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );

  return (*m_coordinates[ dim ])( nodeIdx );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::reserve( IndexType capacity )
{

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
    m_coordinates[ dim ]->reserve( capacity );
  }

  SLIC_ASSERT( consistencyCheck() );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::resize( IndexType size )
{

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
    m_coordinates[ dim ]->resize( size );
  }

  SLIC_ASSERT( consistencyCheck() );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setResizeRatio( double ratio )
{

  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
    m_coordinates[ dim ]->setResizeRatio( ratio );
  }

  SLIC_ASSERT( consistencyCheck() );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::shrink( )
{

  for ( int dim = 0; dim < m_ndims ; ++dim )
  {
    SLIC_ASSERT( m_coordinates[ dim ] != AXOM_NULLPTR );
    m_coordinates[ dim ]->shrink();
  }

  SLIC_ASSERT( consistencyCheck() );
}

//------------------------------------------------------------------------------
// in-line  private method implementations
//------------------------------------------------------------------------------
inline void MeshCoordinates::initialize( IndexType numNodes,
                                         IndexType maxCapacity )
{
  SLIC_ASSERT( numNodes <= maxCapacity );
  SLIC_ASSERT( !this->invalidDimension() );

  for ( int i=0; i < m_ndims; ++i )
  {
    m_coordinates[ i ] = new Array< double >( numNodes, 1, maxCapacity );
  }

  SLIC_ASSERT( consistencyCheck() );
}


} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESHCOORDINATES_HPP_ */
