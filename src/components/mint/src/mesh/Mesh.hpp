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

#ifndef MINT_MESH_HPP_
#define MINT_MESH_HPP_

#include "axom/Macros.hpp"            // for Axom macros

#include "mint/config.hpp"            // for mint compile-time type definitions
#include "mint/FieldAssociation.hpp"  // for FieldAssociation enum
#include "mint/FieldData.hpp"         // for mint::FieldData
#include "mint/MeshType.hpp"          // for MeshType enum and property traits

#include "slic/slic.hpp"              // for SLIC macros


namespace axom
{

/* Forward declarations */
#ifdef MINT_USE_SIDRE
namespace sidre
{
class Group;
}
#endif

namespace mint
{

/* Forward declarations */
class Field;

class Mesh
{
public:
  /*!
   * \brief Constructor.
   * \param [in] ndims the number of dimensions
   * \param [in] type the mesh type.
   * \param [in] blockId the block ID for this mesh instance.
   * \param [in] partId the partition ID for this mesh instance.
   */
  Mesh( int ndims, int type, int blockId, int partId );

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Constructor for use with a group that already has data.
   * \param [in] group the sidre::Group to use.
   * \pre group != AXOM_NULLPTR.
   */
  Mesh( sidre::Group* group );

  /*!
   * \brief Constructor for use with an empty group.
   * \param [in] group the sidre::Group to use.
   * \param [in] ndims the number of dimensions
   * \param [in] type the mesh type.
   * \param [in] blockId the block ID for this mesh instance.
   * \param [in] partId the partition ID for this mesh instance.
   * \pre group != AXOM_NULLPTR.
   */
  Mesh( sidre::Group* group, int ndims, int type, int blockId, int partId );
#endif

  /*!
   * \brief Destructor.
   */
  virtual ~Mesh();

  /*!
   * \brief Returns the dimension for this mesh instance.
   * \return ndims the dimension of this mesh instance.
   * \post ndims >= 1 && ndims <= 3
   */
  inline int getDimension() const
  { return m_ndims; }

  /*!
   * \brief Returns the ID of this mesh instance.
   * \return Id the ID of the mesh.
   */
  inline int getBlockId() const
  { return m_block_idx; }

  /*!
   * \brief Returns the partition ID of this mesh instance.
   * \return partitionId the partition ID of the mesh.
   */
  inline int getPartitionId() const
  { return m_part_idx; }

  /*!
   * \brief Returns the mesh type of this mesh instance.
   * \return meshType the mesh type
   * \see MeshType
   */
  inline int getMeshType() const
  { return m_type; }

  /*!
   * \brief Checks if this mesh instance has explicit coordinates.
   * \return status true iff the mesh defines coordinates explicitly.
   */
  inline bool hasExplicitCoordinates() const
  { return mesh_properties::explicit_coordinates[ m_type ]; }

  /*!
   * \brief Checks if this mesh instance has explicit connectivity.
   * \return status true iff the mesh defines cell connectivity explicitly.
   */
  inline bool hasExplicitConnectivity() const
  { return mesh_properties::explicit_connectivity[ m_type ]; }

  /*!
   * \brief Checks if the mesh has mixed cell types, e.g., consisting of both
   *  triangle and quad elements or hex,pyramid,prisms and tets in 3-D.
   * \return status true iff the mesh has mixed cell types.
   */
  inline bool hasMixedCellTypes() const
  { return m_type==MINT_UNSTRUCTURED_MIXED_ELEMENT_MESH; }

  /*!
   * \brief Checks if this Mesh instance is associated with a Sidre Group.
   * \return status true if the Mesh is associated with a group in a Sidre
   *  hierarchy, else, false.
   */
  inline bool hasSidreGroup( ) const;

/// \name Methods to Create, Access & Remove Fields from a Mesh
/// @{

  /*!
   * \brief Returns const pointer to the FieldData instance with the specified
   *  mesh field association, e.g., NODE_CENTERED, CELL_CENTERED, etc.
   *
   * \param [in] association the specified mesh field association
   * \return fd pointer to the requested FieldData instance
   *
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   * \post fd != AXOM_NULLPTR
   *
   * \see FieldAssociation
   * \see FieldData
   */
  inline const FieldData* getFieldData( int association ) const;

  /*!
   * \brief Check if a field with the given name and association exists.
   *
   * \param [in] name the name of the field in query.
   * \param [in] association the field association, e.g., NODE_CENTERED, etc.
   *
   * \return status true if the field exists, else, false.
   *
   * \pre name.empty()==false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  inline bool hasField( const std::string& name, int association ) const;

  /*!
   * \brief Creates a new field with the given name and specified mesh field
   *  association, e.g., NODE_CENTERED, CELL_CENTERED, etc.
   *
   * \param [in] name the name of the new field.
   * \param [in] association the mesh field association.
   * \param [in] num_components number of components of the field (optional).
   * \param [in] storeInSidre indicates whether to store the field in the
   *  corresponding Sidre group (optional).
   * \param [in] capacity
   *
   * \return ptr raw pointer to the data buffer of the new field.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   * \post ptr != AXOM_NULLPTR
   * \post hasField( name ) == true
   *
   * \see FieldAssociation
   */
  template < typename T >
  inline T* createField( const std::string& name,
                         int association,
                         IndexType num_components=1,
                         bool storeInSidre=true,
                         IndexType capacity=USE_DEFAULT );

  /*!
   * \brief Creates a new field from an external buffer that has the given name
   *  and specified mesh field association, e.g., NODE_CENTERED, CELL_CENTERED,
   *  etc.
   *
   * \param [in] name the name of the new field.
   * \param [in] association the mesh field association.
   * \param [in] data pointer to the external data buffer.
   * \param [in] num_components number of components of the field (optional).
   *
   * \return ptr raw pointer to the data buffer of the new field.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == false
   * \pre data != AXOM_NULLPTR
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   * \post ptr != AXOM_NULLPTR
   * \post ptr == data
   * \post hasField( name ) == true
   *
   * \see FieldAssociation
   */
  template < typename T >
  inline T* createField( const std::string& name,
                         int association,
                         T* data,
                         IndexType num_components=1 );

  /*!
   * \brief Removes the field with the given name and specified association.
   *
   * \param [in] name the name of the field to remove.
   * \param [in] association the mesh field association.
   *
   * \return status true if the field is removed successfully, else, false.
   *
   * \pre name.emtpy() == false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  inline bool removeField( const std::string& name, int association );

  /*!
   * \brief Returns pointer to buffer of the field with the given ane and
   *  specified mesh field association.
   *
   * \param [in] name the name of the requested field.
   * \param [in] association the mesh field association.
   * \param [out] num_components the number of components per tuple (optional).
   *
   * \return ptr raw pointer to the data buffer of the requested field.
   *
   * \pre name.empty() == false
   * \pre hasField( name )
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  /// @{

  template < typename T >
  inline T* getFieldPtr( const std::string& name,
                         int association,
                         IndexType& num_components );

  template < typename T >
  inline T* getFieldPtr( const std::string& name,
                         int association );

  template < typename T >
  inline const T* getFieldPtr( const std::string& name,
                               int association,
                               IndexType& num_components ) const;

  template < typename T >
  inline const T* getFieldPtr( const std::string& name,
                               int association ) const;

  /// @}

/// @}

  /// \name Virtual API
  /// @{

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual IndexType getMeshNumberOfNodes() const = 0;

  /*!
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual IndexType getMeshNumberOfCells() const = 0;

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual int getMeshNumberOfCellNodes( IndexType cellIdx ) const = 0;

  /*!
   * \brief Returns the cell connectivity of the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \param [out] cell user-supplied buffer to store cell connectivity info.
   * \note cell must have sufficient size to hold the connectivity information.
   * \pre cellIdx >= 0 && cellIdx < this->getMeshNumberOfCells()
   * \pre cell != AXOM_NULLPTR.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshCell( IndexType cellIdx, IndexType* cell ) const = 0;

  /*!
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   */
  virtual int getMeshCellType( IndexType cellIdx ) const = 0;

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] coordinates user-supplied buffer to store the node coordinates.
   * \pre 0 <= nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshNode( IndexType nodeIdx,
                            double* coordinates ) const = 0;

  /*!
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < m_ndims
   */
  virtual double getMeshNodeCoordinate( IndexType nodeIdx, int dim ) const = 0;

  /// @}

protected:

  int m_ndims;          /*! mesh dimension */
  int m_type;           /*! the type of the mesh */
  int m_block_idx;      /*! the Block ID of the mesh */
  int m_part_idx;       /*! the partition ID of the mesh */

private:

  /*!
   * \brief Returns the number of tuples for the given mesh field association.
   * \param [in] association the mesh field association, e.g., NODE_CENTERED
   * \return N the number of tuples
   * \post N >= 0
   */
  inline IndexType getNumTuples( int association ) const;

  /*!
   * \brief Allocates the FieldData internal data-structures.
   * \note Helper method that is called from the constructor.
   */
  void allocateFieldData( );

  /*!
   * \brief Deallocates the FieldData internal data-structures.
   * \note Helper method that is called by the destructor.
   */
  void deallocateFieldData( );

  FieldData* m_mesh_fields[ NUM_FIELD_ASSOCIATIONS ];

#ifdef MINT_USE_SIDRE
  sidre::Group * m_group;
#endif

  IndexType * m_num_cells;       /*! The number of cells in the mesh */
  IndexType* m_cell_capacity;    /*! The cell storage capacity */
  double* m_cell_resize_ratio;    /*! The cell resize ratio */

  IndexType* m_num_faces;        /*! The number of faces in the mesh */
  IndexType* m_face_capacity;    /*! The face storage capacity */
  double* m_face_resize_ratio;    /*! The face resize ratio */

  IndexType* m_num_edges;        /*! The number of edges in the mesh */
  IndexType* m_edge_capacity;    /*! The edge storage capacity */
  double* m_edge_resize_ratio;    /*! The edge resize ratio */

  IndexType* m_num_nodes;        /*! The number of nodes in the mesh */
  IndexType* m_node_capacity;    /*! The node storage capacity */
  double* m_node_resize_ratio;    /*! The node resize ratio */


  DISABLE_COPY_AND_ASSIGNMENT( Mesh );
  DISABLE_MOVE_AND_ASSIGNMENT( Mesh );
};

//------------------------------------------------------------------------------
//  IMPLEMENTATION OF TEMPLATE & IN-LINE METHODS
//------------------------------------------------------------------------------

inline bool Mesh::hasSidreGroup( ) const
{
#ifdef MINT_USE_SIDRE
  return ( m_group != AXOM_NULLPTR );
#else
  return false;
#endif
}

//------------------------------------------------------------------------------
inline IndexType Mesh::getNumTuples( int association ) const
{
  IndexType num_tuples = 0;

  switch ( association )
  {
  case NODE_CENTERED:
    num_tuples = *m_num_nodes;
    break;
  case CELL_CENTERED:
    num_tuples = *m_num_cells;
    break;
  case FACE_CENTERED:
    num_tuples = *m_num_faces;
    break;
  default:
    SLIC_ASSERT( association==EDGE_CENTERED );
    num_tuples = *m_num_edges;
  } // END switch

  return ( num_tuples );
}

//------------------------------------------------------------------------------
inline const FieldData* Mesh::getFieldData( int association ) const
{
  SLIC_ERROR_IF( association < 0 || association >= NUM_FIELD_ASSOCIATIONS,
                  "invalid field association [" << association << "]" );
  SLIC_ERROR_IF( m_mesh_fields[ association ]==AXOM_NULLPTR,
              "null field data object w/association [" << association << "]" );
  return m_mesh_fields[ association ];
}

//------------------------------------------------------------------------------
inline bool Mesh::hasField( const std::string& name, int association ) const
{
  const FieldData* fd = getFieldData( association );
  SLIC_ASSERT( fd != AXOM_NULLPTR );
  return fd->hasField( name );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::createField( const std::string& name,
                             int association,
                             IndexType num_components,
                             bool storeInSidre,
                             IndexType capacity )
{
  FieldData* fd = const_cast< FieldData* >( getFieldData( association ) );
  SLIC_ASSERT( fd != AXOM_NULLPTR );

  IndexType num_tuples = getNumTuples( association );
  T* ptr = fd->createField< T >( name, num_tuples, num_components,
                                 storeInSidre, capacity );
  SLIC_ASSERT( ptr != AXOM_NULLPTR );

  return ( ptr );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::createField( const std::string& name,
                             int association,
                             T* data,
                             IndexType num_components )
{
  SLIC_ASSERT( data != AXOM_NULLPTR );

  FieldData* fd = const_cast< FieldData* >( getFieldData( association ) );
  SLIC_ASSERT( fd != AXOM_NULLPTR );

  IndexType num_tuples = getNumTuples( association );
  T* ptr = fd->createField< T >( name, data, num_tuples, num_components );
  SLIC_ASSERT( ptr != AXOM_NULLPTR );

  return ( ptr );
}

//------------------------------------------------------------------------------
inline bool Mesh::removeField( const std::string& name, int association )
{
  bool status   = false;
  FieldData* fd = const_cast< FieldData* >( getFieldData( association ) );

  const bool hasField = fd->hasField( name );
  SLIC_WARNING_IF( !hasField, "field [" << name << "] does not exist!" );

  if ( hasField )
  {
    fd->removeField( name );
    status = true;
  }

  return ( status );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::getFieldPtr( const std::string& name,int association )
{
  IndexType num_components = 0;
  return getFieldPtr< T >( name, association, num_components );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::getFieldPtr( const std::string& name,
                             int association,
                             IndexType& num_components )
{
  const T* ptr = getFieldPtr< T >( name, association, num_components );
  return ( const_cast< T* >( ptr ) );
}

//------------------------------------------------------------------------------
template < typename T >
inline const T* Mesh::getFieldPtr( const std::string& name,
                                   int association ) const
{
  IndexType num_components = 0;
  return getFieldPtr< T >( name, association, num_components );
}

//------------------------------------------------------------------------------
template < typename T >
inline const T* Mesh::getFieldPtr( const std::string& name,
                             int association,
                             IndexType& num_components ) const
{
  const FieldData* fd = getFieldData( association );
  SLIC_ASSERT( fd != AXOM_NULLPTR );

  IndexType num_tuples = 0;
  const T* ptr = fd->getFieldPtr< T >( name, num_tuples, num_components );
  SLIC_ASSERT( ptr != AXOM_NULLPTR );

  return ( ptr );
}

} /* namespace mint */
} /* namespace axom */

#endif /* MESH_HXX_ */
