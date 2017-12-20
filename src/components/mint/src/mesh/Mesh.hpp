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

#ifndef MESH_HXX_
#define MESH_HXX_

#include "axom/Macros.hpp"
#include "mint/FieldData.hpp"
#include "mint/MeshType.hpp"
#include "mint/DataTypes.hpp"


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
   * \brief Returns the FieldData instance associated with cell-centered fields.
   * \return fd pointer to the FieldData instance for cell-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline FieldData& getCellFieldData()
  { return m_cell_data; }

  /*!
   * \brief Returns the FieldData instance associated with cell-centered fields.
   * \return fd pointer to the FieldData instance for cell-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline const FieldData& getCellFieldData() const
  { return m_cell_data; }

  /*!
   * \brief Returns the FieldData instance associated with face-centered fields.
   * \return fd pointer to the FieldData instance for face-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline FieldData& getFaceFieldData()
  { return m_face_data; }

  /*!
   * \brief Returns the FieldData instance associated with face-centered fields.
   * \return fd pointer to the FieldData instance for face-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline const FieldData& getFaceFieldData() const
  { return m_face_data; }

  /*!
   * \brief Returns the FieldData instance associated with edge-centered fields.
   * \return fd pointer to the FieldData instance for edge-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline FieldData& getEdgeFieldData()
  { return m_edge_data; }

  /*!
   * \brief Returns the FieldData instance associated with edge-centered fields.
   * \return fd pointer to the FieldData instance for edge-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline const FieldData& getEdgeFieldData() const
  { return m_edge_data; }

  /*!
   * \brief Returns the FieldData instance associated with node-centered fields.
   * \return fd pointer to the FieldData instance for node-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline FieldData& getNodeFieldData()
  { return m_node_data; }

  /*!
   * \brief Returns the FieldData instance associated with node-centered fields.
   * \return fd pointer to the FieldData instance for node-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline const FieldData& getNodeFieldData() const
  { return m_node_data; }

  /*!
   * \brief Add a cell centered field to the mesh.
   * \param name the name of the field.
   * \param num_components the number of components per value.
   * \tparam T the type of field to add.
   * \return true iff the field was successfully added.
   */
//  template < typename FieldType >
//  Field * addCellField( const std::string& name, int num_components=1 )
//  {
//    int size = getMeshNumberOfCells();
//    int capacity = getMeshCellCapacity();
//    double resize_ratio = getMeshCellResizeRatio();
//    return m_cell_data.addField< FieldType >( name, size, capacity,
//                                              num_components,
//                                              resize_ratio );
//  }
//
//  /*!
//   * \brief Add a face centered field to the mesh.
//   * \param name the name of the field.
//   * \param num_components the number of components per value.
//   * \tparam T the type of field to add.
//   * \return true iff the field was successfully added.
//   */
//  template < typename FieldType >
//  Field * addFaceField( const std::string& name, int num_components=1 )
//  {
//    int size = getMeshNumberOfFaces();
//    int capacity = getMeshCellCapacity();
//    double resize_ratio = getMeshCellResizeRatio();
//    return m_face_data.addField< FieldType >( name, size, capacity,
//                                              num_components,
//                                              resize_ratio );
//  }
//
//  /*!
//   * \brief Add a edge centered field to the mesh.
//   * \param name the name of the field.
//   * \param num_components the number of components per value.
//   * \tparam T the type of field to add.
//   * \return true iff the field was successfully added.
//   */
//  template < typename FieldType >
//  Field * addEdgeField( const std::string& name, int num_components=1 )
//  {
//    int size = getMeshNumberOfEdges();
//    int capacity = getMeshCellCapacity();
//    double resize_ratio = getMeshCellResizeRatio();
//    return m_edge_data.addField< FieldType >( name, size, capacity,
//                                              num_components,
//                                              resize_ratio );
//  }
//
//  /*!
//   * \brief Add a node centered field to the mesh.
//   * \param name the name of the field.
//   * \param num_components the number of components per value.
//   * \tparam T the type of field to add.
//   * \return true iff the field was successfully added.
//   */
//  template < typename FieldType >
//  Field * addNodeField( const std::string& name, int num_components=1 )
//  {
//    int size = getMeshNumberOfNodes();
//    int capacity = getMeshNodeCapacity();
//    double resize_ratio = getMeshNodeResizeRatio();
//    return m_node_data.addField< FieldType >( name, size, capacity,
//                                              num_components,
//                                              resize_ratio );
//  }

  /// \name Virtual API
  /// @{

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual localIndex getMeshNumberOfNodes() const = 0;

//
//  virtual localIndex getMeshNodeCapacity() const = 0;
//
//
//  virtual double getMeshNodeResizeRatio() const = 0;

  /*!
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual localIndex getMeshNumberOfCells() const = 0;

//
//  virtual localIndex getMeshCellCapacity() const = 0;
//
//
//  virtual double getMeshCellResizeRatio() const = 0;
//
//
//  virtual localIndex getMeshNumberOfFaces() const = 0;
//
//
//  virtual localIndex getMeshNumberOfEdges() const = 0;

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual int getMeshNumberOfCellNodes( localIndex cellIdx ) const = 0;

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
  virtual void getMeshCell( localIndex cellIdx, localIndex* cell ) const = 0;

  /*!
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   */
  virtual int getMeshCellType( localIndex cellIdx ) const = 0;

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] coordinates user-supplied buffer to store the node coordinates.
   * \pre 0 <= nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshNode( localIndex nodeIdx,
                            double* coordinates ) const = 0;

  /*!
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < m_ndims
   */
  virtual double getMeshNodeCoordinate( localIndex nodeIdx, int dim ) const = 0;

  /// @}

protected:

  inline void setCellDataSize( localIndex size )
  { m_cell_data.setSize( size ); }


  inline void setFaceDataSize( localIndex size )
  { m_face_data.setSize( size ); }


  inline void setEdgeDataSize( localIndex size )
  { m_edge_data.setSize( size ); }


  inline void setNodeDataSize( localIndex size )
  { m_node_data.setSize( size ); }


  inline void setCellDataCapacity( localIndex capacity )
  { m_cell_data.setCapacity( capacity ); }


  inline void setFaceDataCapacity( localIndex capacity )
  { m_face_data.setCapacity( capacity ); }


  inline void setEdgeDataCapacity( localIndex capacity )
  { m_edge_data.setCapacity( capacity ); }


  inline void setNodeDataCapacity( localIndex capacity )
  { m_node_data.setCapacity( capacity ); }


  inline void setCellDataResizeRatio( double ratio )
  { m_cell_data.setResizeRatio( ratio ); }


  inline void setFaceDataResizeRatio( double ratio )
  { m_face_data.setResizeRatio( ratio ); }


  inline void setEdgeDataResizeRatio( double ratio )
  { m_edge_data.setResizeRatio( ratio ); }


  inline void setNodeDataResizeRatio( double ratio )
  { m_node_data.setResizeRatio( ratio ); }


  int m_ndims;          /*! mesh dimension */
  int m_type;           /*! the type of the mesh */
  int m_block_idx;      /*! the Block ID of the mesh */
  int m_part_idx;       /*! the partition ID of the mesh */

  FieldData m_cell_data; /*! FieldData instance for cell-centered fields. */
  FieldData m_face_data; /*! FieldData instance for face-centered fields. */
  FieldData m_edge_data; /*! FieldData instance for edge-centered fields. */
  FieldData m_node_data; /*! FieldData instance for node-centered fields. */

private:

#ifdef MINT_USE_SIDRE
  sidre::Group * m_group;
#endif

  localIndex * m_num_cells;       /*! The number of cells in the mesh */
  localIndex* m_cell_capacity;    /*! The cell storage capacity */
  double* m_cell_resize_ratio;    /*! The cell resize ratio */

  localIndex* m_num_faces;        /*! The number of faces in the mesh */
  localIndex* m_face_capacity;    /*! The face storage capacity */
  double* m_face_resize_ratio;    /*! The face resize ratio */

  localIndex* m_num_edges;        /*! The number of edges in the mesh */
  localIndex* m_edge_capacity;    /*! The edge storage capacity */
  double* m_edge_resize_ratio;    /*! The edge resize ratio */

  localIndex* m_num_nodes;        /*! The number of nodes in the mesh */
  localIndex* m_node_capacity;    /*! The node storage capacity */
  double* m_node_resize_ratio;    /*! The node resize ratio */


  DISABLE_COPY_AND_ASSIGNMENT(Mesh);
  DISABLE_MOVE_AND_ASSIGNMENT(Mesh);
};

} /* namespace mint */
} /* namespace axom */

#endif /* MESH_HXX_ */
