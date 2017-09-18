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

#ifndef STRUCTUREDMESH_HXX_
#define STRUCTUREDMESH_HXX_

#include "mint/Mesh.hpp"
#include "mint/Extent.hpp"

#include "axom/Types.hpp"
#include "axom/Macros.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <cstddef> // for AXOM_NULLPTR

namespace axom
{
namespace mint
{

class StructuredMesh : public Mesh
{
public:

  /*!
   * \brief Destructor.
   */
  virtual ~StructuredMesh();

  /// \name Virtual Mesh API
  /// @{

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfNodes() const { return this->getNumberOfNodes(); };

  /*!
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfCells() const { return this->getNumberOfCells(); };

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual int getMeshNumberOfCellNodes( int AXOM_NOT_USED( cellIdx ) ) const
  { return this->getNumberOfCellNodes(); }

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
  virtual void getMeshCell( int cellIdx, int* cell ) const
  { this->getCell( cellIdx, cell ); };

  /*!
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   */
  virtual int getMeshCellType( int AXOM_NOT_USED( cellIdx ) ) const;

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [out] coordinates user-supplied buffer to store the node
   *  coordinates.
   * \pre nodeIdx >= && nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshNode( int nodeIdx, double* coordinates ) const
  { this->getNode( nodeIdx, coordinates ); };

  /*!
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < this->getDimension()
   */
  virtual double getMeshNodeCoordinate( int nodeIdx, int dim ) const
  { return (this->getNodeCoordinate( nodeIdx, dim ) ); };

  /// @}

  /// \name Structured Mesh API
  /// @{

  /*!
   * \brief Returns the dimensions of this mesh instance.
   * \param [in] ndims 3-tuple holding the number of nodes in each dimension.
   * \post ndims[ i ] >= 1, \$ \forall i \in [0,2] \$
   */
  inline void getExtentSize( int ndims[ 3 ] ) const
  {
    for ( int i=0 ; i < 3 ; ++i )
    {
      ndims[ i ] = m_extent->size( i );
    } // END for

  }

  /*!
   * \brief Returns stride to the second dimension.
   * \return jp stride to the second dimension.
   * \post jp >= 0.
   */
  inline int jp() const
  { return m_extent->jp(); }

  /*!
   * \brief Returns stride to the third dimension.
   * \return kp stride to the third dimension.
   * \post kp >= 0.
   */
  inline int kp() const
  { return m_extent->kp(); }

  /*!
   * \brief Returns the number of nodes in this mesh instance.
   * \return N the total number of nodes in the mesh.
   * \pre m_extent != AXOM_NULLPTR
   * \post N >= 0.
   */
  inline int getNumberOfNodes() const
  { return m_extent->getNumNodes(); };

  /*!
   * \brief Returns the number of cells in this mesh instance.
   * \return N the total number of cells in the mesh.
   * \pre m_extent != AXOM_NULLPTR.
   * \post N >= 0.
   */
  inline int getNumberOfCells() const
  { return m_extent->getNumCells(); };

  /*!
   * \brief Returns the number of cell nodes.
   * \return N the number of nodes per cell.
   * \post N=4 if 2-D, else N=8.
   */
  int getNumberOfCellNodes() const;

  /*!
   * \brief Returns the linear index corresponding to the given logical indices.
   * \param [in] i logical index of the first dimension.
   * \param [in] j logical index of the second dimension.
   * \param [in] k logical index of the third dimension (optional)
   * \return idx the corresponding linear index.
   */
  inline int getLinearIndex( int i, int j, int k  ) const
  { return m_extent->getLinearIndex( i, j, k ); };

  /*!
   * \brief Returns the linear index corresponding to the given logical indices.
   * \param [in] i logical index of the first dimension.
   * \param [in] j logical index of the second dimension.
   * \return idx the corresponding linear index.
   */
  inline int getLinearIndex( int i, int j ) const
  { return m_extent->getLinearIndex( i, j ); };

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * cell
   * indices.
   * \param [in] i logical cell index of the first dimension.
   * \param [in] j logical cell index of the second dimension.
   * \param [in] k logical cell index of the third dimension (optional)
   * \return idx the corresponding linear index of the cell.
   */
  inline int getCellLinearIndex( int i, int j, int k ) const
  { return m_extent->getCellLinearIndex( i, j, k); };

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * cell
   * indices.
   * \param [in] i logical cell index of the first dimension.
   * \param [in] j logical cell index of the second dimension.
   * \param [in] k logical cell index of the third dimension (optional)
   * \return idx the corresponding linear index of the cell.
   */
  inline int getCellLinearIndex( int i, int j ) const
  { return m_extent->getCellLinearIndex( i, j ); };

  /*!
   * \brief Returns the cell connectivity for the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre cellIdx >= 0 && cellIdx < this->getNumberOfCells()
   * \pre the user-supplied cell buffer must be of getNumberOfCellNodes() size.
   */
  inline void getCell( int cellIdx, int* cell ) const;

  /*!
   * \brief Returns the cell connectivity of the cell at (i,j)
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre this->getDimension() == 2.
   */
  inline void getCell( int i, int j, int* cell ) const;

  /*!
   * \brief Returns the cell connectivity of the cell at (i,j,k)
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [in] k logical index of the cell along the third dimension.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre this->getDimension() == 3.
   */
  inline void getCell( int i, int j, int k, int* cell) const;

  /// \name GetNode() methods -- implemented in concrete instances.
  /// @{

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre coordinates != AXOM_NULLPTR.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes().
   */
  virtual void getNode( int nodeIdx, double* coordinates ) const = 0;

  /*!
   * \brief Returns the coordinates of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 2
   */
  virtual void getNode( int i, int j, double* coordinates ) const = 0;

  /*!
   * \brief Returns the coordinates of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [out] coordinates pointer to buffer to populate with coordinates.
   * \pre this->getDimension() == 3
   */
  virtual void getNode( int i, int j, int k, double* coordinates ) const = 0;

  /*!
   * \brief Returns the coordinate of the given node.
   * \param [in] nodeIdx index of the node in query.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre nodeIdx >= 0 && nodeIdx < this->getNumberOfNodes()
   * \pre idim >= 0 && idim < this->getDimension().
   */
  virtual double getNodeCoordinate( int nodeIdx, int idim  ) const = 0;

  /*!
   * \brief Returns the coordinate value of the node at (i,j)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==2.
   * \pre idim >= 0 && idim < this->getDimension().
   */
  virtual double getNodeCoordinate( int i, int j, int idim ) const = 0;

  /*!
   * \brief Returns the coordinate value of the node at (i,j,k)
   * \param [in] i logical index of the node along the first dimension.
   * \param [in] j logical index of the node along the second dimension.
   * \param [in] k logical index of the node along the third dimension.
   * \param [in] idim requested coordinate dimension.
   * \return x the coordinate value of the node.
   * \pre this->getDimension()==3.
   * \pre idim >= 0 && idim < this->getDimension().
   */
  virtual double getNodeCoordinate( int i, int j, int k, int idim ) const = 0;

  /// @}

  /// @}

protected:

  /*!
   * \brief Default constructor.
   * \note Made private to prevent users from calling it.
   */
  StructuredMesh();

  /*!
   * \brief Constructs a structured mesh instance from the given extent.
   * \param [in] ext the structured mesh extent.
   */
  StructuredMesh( int meshType, int ndims, const int ext[6] );

  /*!
   * \brief Constructs a structured mesh instance from the given extent that is
   *  identified by the given blockId and partitionId pair.
   * \param [in] meshType the structured mesh type.
   * \param [in] ext the structured mesh extent.
   * \param [in] blockId the block ID of the mesh.
   * \param [in] partId the partition ID of the mesh.
   */
  StructuredMesh( int meshType, int ndims, const int ext[6], int blockId,
                  int partId );

  Extent< int >* m_extent;  /*!< grid extent */

private:
  DISABLE_COPY_AND_ASSIGNMENT( StructuredMesh );
  DISABLE_MOVE_AND_ASSIGNMENT( StructuredMesh );
};

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------
namespace axom
{
namespace mint
{

inline void StructuredMesh::getCell(int cellIdx, int* cell) const
{
  SLIC_ASSERT(  cell != AXOM_NULLPTR );
  SLIC_ASSERT(  (cellIdx >= 0) && (cellIdx < this->getNumberOfCells() ) );

  const int* offsets_table = m_extent->getCellOffSets();
  const int num_cell_nodes = this->getNumberOfCellNodes();

  // STEP 0: calculate logical indices of the cell's first corner node.
  int jp_minus_1 = m_extent->jp()-1;
  int ii=0, jj=0, kk=0;

  if  ( this->getDimension() == 1 )
  {

    SLIC_ASSERT( num_cell_nodes==2 );
    ii = cellIdx;

  }
  else if ( this->getDimension() == 2 )
  {

    SLIC_ASSERT( num_cell_nodes==4 );

    ii = cellIdx % jp_minus_1;
    jj = cellIdx / jp_minus_1;

  }
  else
  {

    SLIC_ASSERT(  this->getDimension()==3 );
    SLIC_ASSERT(  num_cell_nodes==8 );

    int kp_minus_1 = m_extent->size(1)-1;
    ii = cellIdx % jp_minus_1;
    jj = cellIdx / jp_minus_1 % kp_minus_1;
    kk = cellIdx / ( jp_minus_1 * kp_minus_1 );

  }

  // STEP 1: calculate linear index of corner node from (ii,jj,kk)
  const int n0  = ii + jj*m_extent->jp() + kk*m_extent->kp();

  // STEP 2: Last, use the offsets table to get the all the cell nodes.
  for ( int i=0 ; i < num_cell_nodes ; ++i )
  {

    cell[ i ] = n0 + offsets_table[ i ];

  } // END for all cell nodes

}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCell(int i, int j, int* cell) const
{
  SLIC_ASSERT(  this->getDimension()==2 );
  SLIC_ASSERT(  cell != AXOM_NULLPTR );

  const int num_cell_nodes = this->getNumberOfCellNodes();
  SLIC_ASSERT(  num_cell_nodes == 4 );

  const int* offsets_table = m_extent->getCellOffSets();

  const int n0 = m_extent->getLinearIndex(i,j);

  for ( int d=0 ; d < num_cell_nodes ; ++d )
  {

    cell[ d ] = n0 + offsets_table[ d ];

  }

}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCell(int i, int j, int k, int* cell) const
{
  SLIC_ASSERT(  this->getDimension()==3 );
  SLIC_ASSERT(  cell != AXOM_NULLPTR );

  const int num_cell_nodes = this->getNumberOfCellNodes();
  SLIC_ASSERT(  num_cell_nodes == 8 );

  const int* offsets_table = m_extent->getCellOffSets();

  const int n0 = m_extent->getLinearIndex(i,j,k);

  for ( int d=0 ; d < num_cell_nodes ; ++d )
  {

    cell[ d ] = n0 + offsets_table[ d ];

  }

}

} /* namespace mint */
} /* namespace axom */

#endif /* STRUCTUREDMESH_HXX_ */
