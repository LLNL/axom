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

#ifndef MINT_UNIFORMMESH_HPP_
#define MINT_UNIFORMMESH_HPP_

#include "mint/config.hpp"                // for IndexType, int64
#include "mint/StructuredMesh.hpp"        // for StructuredMesh base class

#include "slic/slic.hpp"                  // for SLIC macros

namespace axom
{
namespace mint
{

/*!
 * \class UniformMesh
 *
 * \brief Provides the ability to represent and operate on a UniformMesh
 *
 *  The UniformMesh extends from the StructuredMesh base class and provides the
 *  ability to represent and operate on a uniform mesh and associated data.
 *
 *  A <em> uniform mesh </em>, also called a regular mesh, subdivides the domain
 *  in cells that have uniform spacing across each coordinate axis.The nodes and
 *  cells of a <em> uniform </em> mesh are arranged on a regular grid lattice,
 *  where the topology is implicitly defined according to the ordering of the
 *  nodes in the logical i-j-k index space. Moreover, the geometry is also
 *  implicitly defined on a <em> uniform mesh </em>. Given an origin point,
 *  \f$ \hat{x_0} \f$ and spacing, \f$ \hat{h} \f$, along each axis, the
 *  coordinates of a node can be evaluated algebraically by the following:
 *
 *    \f$ \hat{p} = \hat{x_0} + \hat{i} \times \hat{h} \f$,
 *
 *  where \f$\hat{i}\f$ is the logical i-j-k index of the corresponding node.
 *
 *  A UniformMesh object may be constructed using (a) a native constructor, or
 *  (b) from a a Sidre group when Mint is compiled with Sidre support.
 *
 *  * <b> Native Constructor </b> <br />
 *
 *    When using native constructor, the UniformMesh object owns all memory
 *    that is associated with the mesh. Once the UniformMesh object goes
 *    out-of-scope all memory associated with it is returned to the system.
 *
 *  * <b> Sidre Constructor </b> <br />
 *
 *    When a UniformMesh is associated with a Sidre Group, Sidre owns all the
 *    memory. Once the UniformMesh goes out-of-scope the data remains persistent
 *    in Sidre.
 *
 * \note When using Sidre, the specified group must conform to the conventions
 *  defined by the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
 *  computational mesh blueprint </a>
 *
 * \see StructuredMesh
 * \see Extent
 */
class UniformMesh : public StructuredMesh
{
public:

  /*!
   * \brief Default constructor. Disabled.
   */
  UniformMesh() = delete;

/// \name Native Constructors
/// @{

  /*!
   * \brief Constructs a uniform mesh with the given origin, specified spacing
   *  and a mesh extent, which specifies the number of nodes and cells along
   *  each dimension.
   *
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] origin pointer to buffer with the coordinates of the origin
   * \param [in] h pointer to buffer with the spacing on each dimension
   * \param [in] ext the extent of this mesh instance.
   *
   * \note The supplied pointers for `origin` and spacing, `h` must have at
   *  least N entries, where N is the dimension of the uniform mesh.
   *
   *  \note The supplied `ext` pointer must point to a buffer that has at least
   *  \f$ 2 \times N \f$ entries, where N is the dimension of the uniform mesh
   *  given in the following order: [imin, imax, jmin, jmax, kmin, kmax]
   *
   * \pre 1 <= dimension <= 3
   * \pre origin != AXOM_NULLPTR
   * \pre h != AXOM_NULLPTR
   * \pre ext != AXOM_NULLPTR
   *
   * \post getDimesion()==dimension
   * \post getOrigin() != AXOM_NULLPTR
   * \post getSpacing() != AXOM_NULLPTR
   * \post getOrigin()[ i ] == origin[ i ] \f$ \forall i \f$
   * \post getSpacing()[ i ] == h[ i ] \f$ \forall i \f$
   */
  UniformMesh( int dimension,
               const double* origin,
               const double* h,
               const int64* ext );

  /*!
   * \brief Constructs a uniform mesh within a specified rectangular region,
   *  defined by its lower and upper corner points, and a mesh extent which
   *  specifies the number of nodes and cells along each dimension.
   *
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] ext the extent of this mesh instance.
   * \param [in] lower_bound lower corner coordinates of rectangular region.
   * \param [in] upper_bound upper corner coordinates of rectangular region.
   *
   * \note The supplied `lower_bound` and `upper_bound` must point to buffer
   *  that have at least N entries, where N is the dimension of the mesh.
   *
   * \note The supplied `ext` pointer must point to a buffer that has at least
   * \f$ 2 \times N \f$ entries, where N is the dimension of the uniform mesh
   * given in the following order: [imin, imax, jmin, jmax, kmin, kmax]
   *
   * \pre 1 <= dimension <= 3
   * \pre ext != AXOM_NULLPTR
   * \pre lower_bound != AXOM_NULLPTR
   * \pre upper_bound != AXOM_NULLPTR
   *
   * \post getDimesion()==dimension
   * \post getOrigin() != AXOM_NULLPTR
   * \post getSpacing() != AXOM_NULLPTR
   * \post getOrigin()[ i ] == lower_bound[ i ] \f$ \forall i \f$
   */
  UniformMesh( int dimension,
               const int64* ext,
               const double* lower_bound,
               const double* upper_bound  );

  /*!
   * \brief Constructs a uniform mesh within a specified rectangular region,
   *  defined by its lower and upper corner points, and mesh dimensions, Ni,
   *  Nj, Nk.
   *
   * \param [in] dimension the dimension of this mesh instance.
   * \param [in] lower_bound lower corner coordinates of rectangular region.
   * \param [in] upper_bound upper cordner coordinates of rectangular region.
   * \param [in] Ni the number of nodes in the i-direction
   * \param [in] Nj the number of nodes in the j-direction (if dimension >= 2)
   * \param [in] Nk the number of nodes in the k-direction (if dimension == 3)
   *
   * \pre 1 <= dimension <= 3
   * \pre lower_bound != AXOM_NULLPTR
   * \pre upper_bound != AXOM_NULLPTR
   * \pre Ni >= 1
   * \pre Nj >= 1 iff dimension >= 2
   * \pre Nk >= 1 iff dimension == 3
   *
   * \post getDimension() == dimension
   * \post getOrigin() != AXOM_NULLPTR
   * \post getSpacing() != AXOM_NULLPTR
   * \post getOrigin()[ i ] == lower_bound[ i ] \f$ \forall i \f$
   */
  UniformMesh( int dimension,
               const double* lower_bound,
               const double* upper_bound,
               IndexType Ni,
               IndexType Nj=-1,
               IndexType Nk=-1 );
/// @}

#ifdef MINT_USE_SIDRE

/// \name Sidre Constructors
/// @{

  /*!
   * \brief Creates a uniform mesh instance from the given Sidre group that
   *  holds uniform mesh data according to the computational mesh blueprint
   *  conventions.
   *
   * \param [in] group pointer to the root group within a Sidre hierarchy.
   * \param [in] topo the name of the topology for this mesh (optional).
   *
   * \note The supplied group is expected to contain uniform mesh data that is
   *  valid & conforming to the conventions described by conduit's
   *  <a href="http://llnl-conduit.readthedocs.io/en/latest/"> computational
   *  mesh blueprint </a>.
   *
   * \note If a topology name is not provided, the implementation will use the
   *  1st topology group under the parent "topologies" group.
   *
   * \note When using this constructor, all data is owned by Sidre. Once the
   *  mesh object goes out-of-scope, the data will remain persistent in Sidre.
   *
   * \pre group != AXOM_NULLPTR
   * \pre blueprint::isValidRootGroup( group )
   * \post hasSidreGroup() == true
   */
  explicit UniformMesh( sidre::Group* group, const std::string& topo="" );

  /*!
   * \brief Constructs a uniform mesh object, on an empty Sidre group, that
   *  covers a rectangular region, defined by its lower and upper corner points
   *  and the desired mesh extent, which specifies the number of nodes and cells
   *  along each dimension.
   *
   * \param [in] dimension the dimension of the mesh
   * \param [in] lower_bound lower corner coordinates of rectangular region
   * \param [in] upper_bound upper corner coordinates of rectangular region
   * \param [in] extent extent of this mesh instance.
   * \param [in] group pointer to the Sidre group
   * \param [in] topo the name of the topology (optional)
   * \param [in] coordset the name of the corresponding coordset (optional)
   *
   * \note The supplied `lower_bound` and `upper_bound` must point to buffer
   *  that have at least N entries, where N is the dimension of the mesh.
   *
   * \note The supplied `ext` pointer must point to a buffer that has at least
   * \f$ 2 \times N \f$ entries, where N is the dimension of the uniform mesh
   * given in the following order: [imin, imax, jmin, jmax, kmin, kmax]
   *
   * \note If a topology and/or coordset name are not provided, internal
   *  defaults will be used by the implementation.
   *
   * \note When using this constructor, all data is owned by Sidre. Once the
   *  UniformMesh object goes out-of-scope, the data will remain persistent in
   *  Sidre.
   *
   * \pre 1 <= dimension <= 3
   * \pre ext != AXOM_NULLPTR
   * \pre lower_bound != AXOM_NULLPTR
   * \pre upper_bound != AXOM_NULLPTR
   * \pre group != AXOM_NULLPTR
   * \pre group->getNumViews()==0
   * \pre group->getNumGroups()==0
   *
   * \post getDimesion()==dimension
   * \post getOrigin() != AXOM_NULLPTR
   * \post getSpacing() != AXOM_NULLPTR
   * \post getOrigin()[ i ] == lower_bound[ i ] \f$ \forall i \f$
   * \post hasSidreGroup() == true
   */
  UniformMesh( int dimension,
               const double* lower_bound,
               const double* upper_bound,
               const int64* extent,
               sidre::Group* group,
               const std::string& topo="",
               const std::string& coordset="" );

  /*!
   * \brief Constructs a uniform mesh object, on an empty Sidre group, that
   *  covers a rectangular region, defined by its lower and upper corner points,
   *  and desired mesh dimensions, \f$ N_i, N_j, N_k f\$
   *
   * \param [in] dimension the dimension of the mesh
   * \param [in] lower_bound lower corner coordinates of rectangular region
   * \param [in] upper_bound upper corner coordinates of rectangular region
   * \param [in] group pointer to the Sidre group
   * \param [in] topo the name of the topology (optional)
   * \param [in] coordset name of the corresponding coordset group (optional)
   * \param [in] Ni number of nodes in the i-direction
   * \param [in] Nj number of nodes in the j-direction (if dimension >= 2)
   * \param [in] Nk number of nodes in the k-direction (if dimension == 3)
   *
   * \note The supplied `lower_bound` and `upper_bound` must point to buffer
   *  that have at least N entries, where N is the dimension of the mesh.
   *
   * \note When using this constructor, all data is owned by Sidre. Once the
   *  UniformMesh object goes out-of-scope, the data will remain persistent in
   *  Sidre.
   *
   * \note If a topology and/or coordset name are not provided by the caller,
   *  internal defaults will be used by the implementation.
   *
   * \pre 1 <= dimension <= 3
   * \pre lower_bound != AXOM_NULLPTR
   * \pre upper_bound != AXOM_NULLPTR
   * \pre Ni >= 1
   * \pre Nj >= 1 iff dimension >= 2
   * \pre Nk >= 1 iff dimension == 3
   * \pre group != AXOM_NULLPTR
   * \pre group->getNumViews()==0
   * \pre group->getNumGroups()==0
   *
   * \post getDimension() == dimension
   * \post getOrigin() != AXOM_NULLPTR
   * \post getSpacing() != AXOM_NULLPTR
   * \post getOrigin()[ i ] == lower_bound[ i ] \f$ \forall i \f$
   * \post hasSidreGroup() == true
   */
  /// @{
  UniformMesh( int dimension,
               const double* lower_bound,
               const double* upper_bound,
               sidre::Group* group,
               const std::string& topo,
               const std::string& coordset,
               IndexType Ni,
               IndexType Nj=-1,
               IndexType Nk=-1 );

  UniformMesh( int dimension,
               const double* lower_bound,
               const double* upper_bound,
               sidre::Group* group,
               IndexType Ni,
               IndexType Nj=-1,
               IndexType Nk=-1  );
  /// @}

/// @}
#endif

/// \name Virtual methods
/// @{

  /*!
   * \brief Destructor.
   */
  virtual ~UniformMesh() { }

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
  {
    SLIC_ERROR( "getCoordinateArray() is not supported for UniformMesh" );
    return AXOM_NULLPTR;
  }

  virtual
  const double* getCoordinateArray(int AXOM_NOT_USED(dim)) const final override
  {
    SLIC_ERROR("getCoordinateArray() is not supported for UniformMesh" );
    return AXOM_NULLPTR;
  }

  /// @}

/// @}

/// @}

/// \name Attribute Querying Methods
/// @{

  /*!
   * \brief Returns a const pointer to origin of the Uniform Mesh
   * \return origin pointer to the buffer storing the coordinates of the origin
   * \post origin != AXOM_NULLPTR
   */
  const double* getOrigin( ) const
  { return m_origin; }

  /*!
   * \brief Returns a const pointer to spacing of the Uniform Mesh.
   * \return h user-supplied buffer to store the spacing of the mesh.
   */
  const double* getSpacing( ) const
  { return m_h; }

/// @}

  /*!
   * \brief Evaluates the physical coordinate of a node along a given diction.
   *
   * \param [in] i the node's logical grid index along the specified direction
   * \param [in] direction the direction in query, e.g, I_DIRECTION, etc.
   *
   * \return x the physical coordinate of the node in the specified direction.
   *
   * \pre i >= 0 && i < getNumberNodeAlongDim( direction )
   * \pre direction >= 0 && direction < getDimension()
   */
  inline double evaluateCoordinate( IndexType i, int direction ) const
  {
    SLIC_ASSERT( direction >=0 && direction < getDimension() );
    SLIC_ASSERT( (i >= 0) && ( i < getNumberOfNodesAlongDim( direction ) ) );

    const double* x0 = getOrigin( );
    const double* h  = getSpacing( );
    return ( x0[ direction ] + i*h[ direction ] );
  }

private:

  double m_origin[ 3 ] = { 0.0, 0.0, 0.0 };
  double m_h[ 3 ]      = { 1.0, 1.0, 1.0 };

  DISABLE_COPY_AND_ASSIGNMENT(UniformMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(UniformMesh);
};


} /* namespace mint */
} /* namespace axom */

#endif /* UNIFORMMESH_HXX_ */
