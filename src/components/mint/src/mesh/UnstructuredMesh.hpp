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

#ifndef UNSTRUCTUREDMESH_HXX_
#define UNSTRUCTUREDMESH_HXX_

#include "mint/Mesh.hpp"
#include "mint/MeshType.hpp"
#include "mint/CellType.hpp"
#include "mint/CellConnectivity.hpp"
#include "mint/MeshCoordinates.hpp"
#include "mint/Field.hpp"
#include "mint/FieldTypes.hpp"
#include "mint/FieldData.hpp"

#include "slic/slic.hpp"

#include "axom/Macros.hpp"
#include "axom/Types.hpp"

#include <cstring> // for memcpy()
#include <fstream> // for fstream

namespace axom
{
namespace mint
{

template < int CellType >
class UnstructuredMesh : public Mesh
{
public:

  /*!
   * \brief Custom constructor. Creates an unstructured mesh of given dimension.
   * \param [in] ndims the number of dimension.
   */
  explicit UnstructuredMesh(int ndims);

  /*!
   * \brief Creates an unstructured mesh of given dimension that is identified
   *  by the supplied blockId and partitionId pair.
   * \param [in] ndims the number of dimension.
   * \param [in] blockId the block ID of the mesh
   * \param [in] partId the partition ID of the mesh
   */
  UnstructuredMesh(int ndims, int blockId, int partId);

  /*!
   * \brief Destructor.
   */
  virtual ~UnstructuredMesh();

  /// \name Virtual API
  /// @{

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfNodes() const
  { return this->getNumberOfNodes(); };

  /*!
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual int getMeshNumberOfCells() const
  { return this->getNumberOfCells(); };

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual int getMeshNumberOfCellNodes( int cellIdx ) const
  { return this->getNumberOfCellNodes( cellIdx ); };

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
  {
    const int numNodes = this->getNumberOfCellNodes( cellIdx );
    const int* myCell  = this->getCell( cellIdx );
    memcpy( cell, myCell, numNodes*sizeof(int) );
  }

  /*!
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   */
  virtual int getMeshCellType( int cellIdx ) const
  { return m_cell_connectivity->getCellType( cellIdx ); }

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] coordinates user-supplied buffer to store the node coordinates.
   * \pre nodeIdx >= && nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshNode( int nodeIdx, double* coordinates ) const
  {
    for ( int i=0 ; i < m_ndims ; ++i )
    {
      const double* coord = this->getMeshCoordinateArray( i );
      coordinates[ i ] = coord[ nodeIdx ];
    } // END for all dimensions
  }

  /*!
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < m_ndims
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual double getMeshNodeCoordinate( int nodeIdx, int dim ) const
  {
    SLIC_ASSERT( dim >= 0 && dim < m_ndims );
    return this->getMeshCoordinateArray( dim )[ nodeIdx ];
  }

  /// @}

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return N the total number of nodes in the mesh.
   * \pre m_node_coordinates != AXOM_NULLPTR.
   * \post N >= 0.
   */
  int getNumberOfNodes() const
  { return m_node_coordinates->getNumberOfPoints(); };

  /*!
   * \brief Returns the total number cells in the mesh.
   * \return N the total number of cells in the mesh.
   * \pre m_cell_connectivity != AXOM_NULLPTR.
   * \post N >= 0.
   */
  int getNumberOfCells() const
  { return m_cell_connectivity->getNumberOfCells(); };

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return nnodes the number of nodes for the given cell.
   * \pre nnodes >= 1.
   */
  int getNumberOfCellNodes( int cellIdx ) const
  { return m_cell_connectivity->getNumberOfNodes( cellIdx ); };

  /*!
   * \brief Inserts a new cell in the mesh
   * \param [in] cell pointer to the connectivity of the cell.
   * \param [in] cell_type the type of cell to insert.
   * \param [in] numNodes number of nodes for the cell.
   * \pre cell != AXOM_NULLPTR.
   */
  void insertCell( const int* cell, int cell_type, int numNodes );

  /*!
   * \brief Inserts a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \pre m_ndims == 1.
   */
  void insertNode( double x );

  /*!
   * \brief Inserts a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \param [in] y the y--coordinate of the new mesh node.
   * \pre m_ndims == 2.
   */
  void insertNode( double x, double y );

  /*!
   * \brief Inserts a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \param [in] y the y--coordinate of the new mesh node.
   * \param [in] z the z--coordinate of the new mesh node.
   * \pre m_ndims == 3.
   */
  void insertNode( double x, double y, double z );

  /*!
   * \brief Inserts a new node in the mesh.
   * \param [in] node pointer to buffer consisting of the coordinates
   * \pre node != AXOM_NULLPTR
   * \pre node must at least be m_ndims long
   */
  void insertNode( const double* node );

  /*!
   * \brief Returns pointer to coordinates array for the requested dimension.
   * \param [in] dim the requested dimension.
   * \return ptr pointer to the coordinates array for the given dimension.
   * \pre m_coordinates != AXOM_NULLPTR.
   * \pre dim < m_ndims
   * \post ptr != AXOM_NULLPTR.
   */
  const double* getMeshCoordinateArray( int dim ) const;

  /*!
   * \brief Returns pointer to the connectivity array of the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \return cell_ptr pointer to the connectivity array of the cell.
   * \pre m_cell_connectivity != AXOM_NULLPTR
   * \pre cellIdx >= 0 && cellIdx < this->getNumberOfCells()
   * \post cell_ptr != AXOM_NULLPTR.
   */
  const int* getCell( int cellIdx ) const;

private:

  /*!
   * \brief Default constructor.
   * \note Made private to prevent users from calling it.
   */
  UnstructuredMesh();

  MeshCoordinates* m_node_coordinates;
  CellConnectivity< int, CellType >* m_cell_connectivity;

  DISABLE_COPY_AND_ASSIGNMENT(UnstructuredMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(UnstructuredMesh);
};

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
//      UnstructuredMesh Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace mint
{

template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh() :
  Mesh(-1,-1,-1,-1),
  m_node_coordinates( AXOM_NULLPTR ),
  m_cell_connectivity( AXOM_NULLPTR )
{}

//------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh(int ndims ) :
  Mesh(ndims, mesh_properties::mesh_of_cell_type[ CellType ], 0, 0 ),
  m_node_coordinates( new MeshCoordinates( ndims ) ),
  m_cell_connectivity( new CellConnectivity< int, CellType >() )
{}

//------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh(
  int ndims, int blockId, int partId ) :
  Mesh(ndims, mesh_properties::mesh_of_cell_type[CellType],blockId,partId ),
  m_node_coordinates( new MeshCoordinates( ndims ) ),
  m_cell_connectivity( new CellConnectivity< int, CellType >() )
{}

//------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::~UnstructuredMesh()
{
  delete m_node_coordinates;
  delete m_cell_connectivity;
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::insertCell( const int* cell, int cell_type,
                                               int num_nodes )
{
  SLIC_ASSERT(  cell != AXOM_NULLPTR );
  SLIC_ASSERT(  m_cell_connectivity != AXOM_NULLPTR );
  SLIC_ASSERT(  cell::num_nodes[ cell_type ] == num_nodes );
  SLIC_ASSERT(  (CellType == MINT_MIXED_CELL) ? true : CellType == cell_type );

  m_cell_connectivity->insertCell( cell, cell_type, num_nodes );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::insertNode( double x )
{
  SLIC_ASSERT(  m_ndims == 1 );
  SLIC_ASSERT(  m_node_coordinates != AXOM_NULLPTR );
  m_node_coordinates->insertPoint( x );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::insertNode( double x, double y )
{
  SLIC_ASSERT(  m_ndims == 2 );
  SLIC_ASSERT(  m_node_coordinates != AXOM_NULLPTR );
  m_node_coordinates->insertPoint( x, y );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::insertNode( double x, double y, double z )
{
  SLIC_ASSERT(  m_ndims == 3 );
  SLIC_ASSERT(  m_node_coordinates != AXOM_NULLPTR );
  m_node_coordinates->insertPoint( x, y, z );
}

//------------------------------------------------------------------------------
template < int CellType >
inline void UnstructuredMesh< CellType >::insertNode( const double* node )
{
  SLIC_ASSERT(  node != AXOM_NULLPTR );
  SLIC_ASSERT(  m_node_coordinates != AXOM_NULLPTR );

  if ( m_ndims == 2 )
  {

    m_node_coordinates->insertPoint( node[0], node[1] );

  }
  else
  {

    SLIC_ASSERT( m_ndims == 3 );
    m_node_coordinates->insertPoint( node[0], node[1], node[2] );

  }

}

//------------------------------------------------------------------------------
template < int CellType >
inline
const double* UnstructuredMesh< CellType >::getMeshCoordinateArray(int dim)
const
{
  SLIC_ASSERT(  m_node_coordinates != AXOM_NULLPTR );
  SLIC_ASSERT(  dim < m_ndims );
  return m_node_coordinates->getCoordinateArray( dim );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
const int* UnstructuredMesh< CellType >::getCell( int cellIdx ) const
{
  SLIC_ASSERT(  m_cell_connectivity != AXOM_NULLPTR );
  SLIC_ASSERT(  cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
  return (*m_cell_connectivity)[ cellIdx ];
}

} /* namespace mint */
} /* namespace axom */

#endif /* UNSTRUCTUREDMESH_HXX_ */
