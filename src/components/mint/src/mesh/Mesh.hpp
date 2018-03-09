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
#include "mint/MeshType.hpp"

namespace axom
{
namespace mint
{

// Forward Declarations
class FieldData;

class Mesh
{
public:

  /*!
   * \brief Destructor.
   */
  virtual ~Mesh();

  /*!
   * \brief Returns the dimension for this mesh instance.
   * \return ndims the dimension of this mesh instance.
   * \post ndims >= 1 && ndims <= 3
   */
  inline int getDimension() const { return m_ndims; };

  /*!
   * \brief Returns the ID of this mesh instance.
   * \return Id the ID of the mesh.
   */
  inline int getBlockId() const { return m_block_idx; };

  /*!
   * \brief Returns the partition ID of this mesh instance.
   * \return partitionId the partition ID of the mesh.
   */
  inline int getPartitionId() const { return m_part_idx; };

  /*!
   * \brief Returns the mesh type of this mesh instance.
   * \return meshType the mesh type
   * \see MeshType
   */
  inline int getMeshType() const { return m_type; };

  /*!
   * \brief Checks if this mesh instance has explicit coordinates.
   * \return status true iff the mesh defines coordinates explicitly.
   */
  inline bool hasExplicitCoordinates() const
  { return mesh_properties::explicit_coordinates[ m_type ]; };

  /*!
   * \brief Checks if this mesh instance has explicit connectivity.
   * \return status true iff the mesh defines cell connectivity explicitly.
   */
  inline bool hasExplicitConnectivity() const
  { return mesh_properties::explicit_connectivity[ m_type ]; };

  /*!
   * \brief Checks if the mesh has mixed cell types, e.g., consisting of both
   *  triangle and quad elements or hex,pyramid,prisms and tets in 3-D.
   * \return status true iff the mesh has mixed cell types.
   */
  inline bool hasMixedCellTypes() const
  { return m_type==MINT_UNSTRUCTURED_MIXED_ELEMENT_MESH; };

  /*!
   * \brief Returns the FieldData instance associated with node-centered fields.
   * \return fd pointer to the FieldData instance for node-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline FieldData* getNodeFieldData() const { return m_node_data; };

  /*!
   * \brief Returns the FieldData instance associated with cell-centered fields.
   * \return fd pointer to the FieldData instance for cell-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline FieldData* getCellFieldData() const { return m_cell_data; };

  /*!
   * \brief Returns the FieldData instance associated with face-centered fields.
   * \return fd pointer to the FieldData instance for face-centered fields.
   * \post fd != AXOM_NULLPTR
   */
  inline FieldData* getFaceFieldData() const { return m_face_data; };

  /// \name Virtual API
  /// @{

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfNodes() const = 0;

  /*!
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfCells() const = 0;

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual int getMeshNumberOfCellNodes( int cellIdx ) const = 0;

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
  virtual void getMeshCell( int cellIdx, int* cell ) const = 0;

  /*!
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   */
  virtual int getMeshCellType( int cellIdx ) const = 0;

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] coordinates user-supplied buffer to store the node coordinates.
   * \pre 0 <= nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshNode( int nodeIdx, double* coordinates ) const = 0;

  /*!
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < m_ndims
   */
  virtual double getMeshNodeCoordinate( int nodeIdx, int dim ) const = 0;

  /// @}

protected:

  /*!
   * \brief Default Constructor.
   * \note Made protected since this is an abstract class.
   */
  Mesh();

  /*!
   * \brief Custom constructor.
   * \param [in] ndims the number of dimensions
   * \param [in] type the mesh type.
   * \param [in] blockId the block ID for this mesh instance.
   * \param [in] partId the partition ID for this mesh instance.
   */
  Mesh( int ndims, int type, int blockId, int partId );

  int m_ndims;          /*!< mesh dimension */
  int m_type;           /*!< the type of the mesh */
  int m_block_idx;      /*!< the Block ID of the mesh */
  int m_part_idx;       /*!< the partition ID of the mesh */

  FieldData* m_node_data;  /*!< FieldData instance for node-centered fields. */
  FieldData* m_cell_data;  /*!< FieldData instance for cell-centered fields. */
  FieldData* m_face_data;  /*!< FieldData instance for face-centered fields. */

private:
  DISABLE_COPY_AND_ASSIGNMENT(Mesh);
  DISABLE_MOVE_AND_ASSIGNMENT(Mesh);
};

} /* namespace mint */
} /* namespace axom */

#endif /* MESH_HXX_ */
