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
#include "mint/Vector.hpp"
#include "mint/DataTypes.hpp"

#include "slic/slic.hpp"

#include "axom/Macros.hpp"
#include "axom/Types.hpp"

#include <cstring> // for std::memcpy
#include <cmath>   // for std::pow

namespace axom
{
namespace mint
{

template < int CellType >
class UnstructuredMesh : public Mesh
{
public:

//  /*!
//   * \brief Custom constructor. Creates an unstructured mesh of given
// dimension.
//   * \param [in] ndims the number of dimension.
//   * \param [in] nodeCapacity size of node buffer to allocate.
//   */
//  UnstructuredMesh( int ndims, localIndex nodeCapacity=100,
//                             localIndex cellNodeCapacity=0,
//                             double nodeResizeRatio=2.0,
//                             double cellResizeRatio=2.0 );
//
//  /*!
//   * \brief Creates an unstructured mesh of given dimension that is identified
//   *  by the supplied blockId and partitionId pair.
//   * \param [in] ndims the number of dimension.
//   * \param [in] blockId the block ID of the mesh.
//   * \param [in] partId the partition ID of the mesh.
//   * \param [in] nodeCapacity size of node buffer to allocate.
//   */
//  UnstructuredMesh( int ndims, int blockId, int partId,
//                    localIndex nodeCapacity=100, localIndex
// cellNodeCapacity=0,
//                    double nodeResizeRatio=2.0, double cellResizeRatio=2.0 );

  UnstructuredMesh( int ndims );

  UnstructuredMesh( int ndims, localIndex num_nodes, localIndex num_cells );

  /*!
   * \brief Destructor.
   */
  virtual ~UnstructuredMesh()
  {}

  /// \name Virtual API
  /// @{

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual localIndex getMeshNumberOfNodes() const override
  { return getNumberOfNodes(); };


//  virtual localIndex getMeshNodeCapacity() const override
//  { return getNodeCapacity(); }
//
//
//  virtual double getMeshNodeResizeRatio() const override
//  { return getNodeResizeRatio(); }

  /*!
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   */
  virtual localIndex getMeshNumberOfCells() const override
  { return getNumberOfCells(); };


//  virtual localIndex getMeshCellCapacity() const override
//  {
//    double capacity_ratio = getCellNodeCapacity();
//    capacity_ratio /= getNumberOfCells();
//    return std::ceil( capacity_ratio * getNumberOfCells() );
//  }
//
//
//  virtual double getMeshCellResizeRatio() const override
//  { return getCellResizeRatio(); }
//
//
//  virtual localIndex getMeshNumberOfFaces() const override
//  { return 0; }
//
//
//  virtual localIndex getMeshNumberOfEdges() const override
//  { return 0; }


  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual int getMeshNumberOfCellNodes( localIndex cellIdx ) const override
  { return getNumberOfCellNodes( cellIdx ); }

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
  virtual void getMeshCell( localIndex cellIdx,
                            localIndex* cell ) const override
  {
    const int numNodes = getNumberOfCellNodes( cellIdx );
    const localIndex* myCell  = getCell( cellIdx );
    std::memcpy( cell, myCell, numNodes * sizeof(localIndex) );
  }

  /*!
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   */
  virtual int getMeshCellType( localIndex cellIdx ) const override
  { return m_cell_connectivity->getCellType( cellIdx ); }

  /*!
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] coordinates user-supplied buffer to store the node coordinates.
   * \pre nodeIdx >= && nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual void getMeshNode( localIndex nodeIdx,
                            double* coordinates ) const override
  {
    for ( int i=0 ; i < getDimension() ; ++i )
    {
      const double* coord = getMeshCoordinateArray( i );
      coordinates[ i ] = coord[ nodeIdx ];
    } // END for all dimensions
  }

  /*!
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < this->getDimension()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   */
  virtual double getMeshNodeCoordinate( localIndex nodeIdx,
                                        int dim ) const override
  {
    SLIC_ASSERT( dim >= 0 && dim < this->getDimension() );
    return getMeshCoordinateArray( dim )[ nodeIdx ];
  }

  /// @}

  /*!
   * \brief Returns the total number of nodes in the mesh.
   * \return N the total number of nodes in the mesh.
   * \post N >= 0.
   */
  localIndex getNumberOfNodes() const
  { return m_node_coordinates->getSize(); };

  /*!
   * \brief Get the number of nodes that can be stored.
   * \return The capacity of the node array.
   */
//  localIndex getNodeCapacity() const
//  { return m_node_coordinates.getCapacity(); };

//  /*!
//   * \brief Set the number of nodes that can be stored.
//   */
//  void setNodeCapacity( localIndex capacity )
//  {
//    m_node_coordinates.setCapacity( capacity );
//    this->setNodeDataCapacity( capacity );
//  }


//  double getNodeResizeRatio() const
//  { return m_node_coordinates.getResizeRatio(); }


//  void setNodeResizeRatio( double ratio )
//  {
//    m_node_coordinates.setResizeRatio( ratio );
//    this->setNodeDataResizeRatio( ratio );
//  }

  /*!
   * \brief Returns the total number cells in the mesh.
   * \return N the total number of cells in the mesh.
   * \post N >= 0.
   */
  localIndex getNumberOfCells() const
  { return m_cell_connectivity->getNumberOfCells(); };

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return nnodes the number of nodes for the given cell.
   * \pre nnodes >= 1.
   */
  int getNumberOfCellNodes( localIndex cellIdx ) const
  { return m_cell_connectivity->getNumberOfNodes( cellIdx ); }

//  /*!
//   * \brief Get the number of nodes that can be stored in the cell
// connectivity.
//   * \return the number of nodes that can be stored in the cell connectivity.
//   */
//  localIndex getCellNodeCapacity() const
//  { return m_cell_connectivity.getCapacity(); }
//
//  /*!
//   * \brief Set the number of nodes that can be stored in the cell
// connectivity.
//   * \param capacity the number of nodes that can be stored in the cell
//   *  connectivity.
//   */
//  void setCellNodeCapacity( localIndex capacity )
//  { m_cell_connectivity.setCellCapacity( capacity ); }


//  double getCellResizeRatio() const
//  { return m_cell_connectivity.getResizeRatio(); }
//
//
//  void setCellResizeRatio( double ratio )
//  {
//    m_cell_connectivity.setResizeRatio( ratio );
//    this->setCellDataResizeRatio( ratio );
//    this->setFaceDataResizeRatio( ratio );
//    this->setEdgeDataResizeRatio( ratio );
//  }

  /*!
   * \brief Adds a new cell in the mesh
   * \param [in] cell pointer to the connectivity of the cell.
   * \param [in] cell_type the type of cell to Add.
   * \pre cell != AXOM_NULLPTR.
   */
  void addCell( const localIndex* cell, int cell_type );

  /*!
   * \brief Adds a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \pre this->getDimension() == 1.
   */
  void addNode( double x );

  /*!
   * \brief Adds a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \param [in] y the y--coordinate of the new mesh node.
   * \pre this->getDimension() == 2.
   */
  void addNode( double x, double y );

  /*!
   * \brief Adds a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \param [in] y the y--coordinate of the new mesh node.
   * \param [in] z the z--coordinate of the new mesh node.
   * \pre this->getDimension() == 3.
   */
  void addNode( double x, double y, double z );

  /*!
   * \brief Adds a new node in the mesh.
   * \param [in] node pointer to buffer consisting of the coordinates
   * \pre node != AXOM_NULLPTR
   * \pre node must at least be this->getDimension() long
   */
  void addNode( const double* node );

  /*!
   * \brief Returns pointer to coordinates array for the requested dimension.
   * \param [in] dim the requested dimension.
   * \return ptr pointer to the coordinates array for the given dimension.
   * \pre dim < this->getDimension()
   * \post ptr != AXOM_NULLPTR.
   */
  const double* getMeshCoordinateArray( int dim ) const;

  /*!
   * \brief Returns pointer to the connectivity array of the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \return cell_ptr pointer to the connectivity array of the cell.
   * \pre cellIdx >= 0 && cellIdx < this->getNumberOfCells()
   * \post cell_ptr != AXOM_NULLPTR.
   */
  const localIndex* getCell( localIndex cellIdx ) const;

private:

  /*!
   * \brief A helper function used in the constructor. Given the dimension and
   *  a node capacity approximates the capacity of m_cell_connectivity.
   * \param [in] ndims the number of dimensions of the mesh.
   * \param [in] nodeCapacity the number of nodes allocated for the mesh.
   * \param [in] cellNodeCapacity the number of nodes to allocate for
   * m_cell_connectivity.
   * \return ptr the number of nodes to allocate for m_cell_connectivity.
   */
  static localIndex approxNumCellNodes( int ndims, localIndex nodeCapacity,
                                        localIndex cellNodeCapacity )
  {
    if ( cellNodeCapacity > 0 )
    {
      return cellNodeCapacity;
    }

    if ( nodeCapacity <= 1)
    {
      return 0;
    }

    double linear_num_nodes = std::pow( nodeCapacity, 1.0 / ndims );
    double linear_num_cells = linear_num_nodes - 1;
    if ( linear_num_nodes < 0 )
    {
      linear_num_nodes = 0;
    }
    return std::pow( linear_num_cells, ndims ) * cell::num_nodes[ CellType ];
  }

  /*!
   * \brief Default constructor.
   * \note Made private to prevent users from calling it.
   */
  UnstructuredMesh()
  {}

  MeshCoordinates* m_node_coordinates;
  CellConnectivity< CellType >* m_cell_connectivity;

  DISABLE_COPY_AND_ASSIGNMENT( UnstructuredMesh );
  DISABLE_MOVE_AND_ASSIGNMENT( UnstructuredMesh );
};


//------------------------------------------------------------------------------
//      UnstructuredMesh Implementation
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//template < int CellType >
//UnstructuredMesh< CellType >::UnstructuredMesh( int ndims,
//                                                localIndex nodeCapacity,
//                                                localIndex cellNodeCapacity,
//                                                double nodeResizeRatio,
//                                                double cellResizeRatio ) :
//  Mesh( ndims, mesh_properties::mesh_of_cell_type[ CellType ], 0, 0 )
//// TODO: ???
////  m_node_coordinates( ndims, nodeCapacity, nodeResizeRatio ),
////  m_cell_connectivity( approxNumCellNodes( ndims, nodeCapacity,
////                                           cellNodeCapacity ),
////                       cellResizeRatio )
//{}

//------------------------------------------------------------------------------
//template < int CellType >
//UnstructuredMesh< CellType >::UnstructuredMesh( int ndims, int blockId,
//                                                int partId,
//                                                localIndex nodeCapacity,
//                                                localIndex cellNodeCapacity,
//                                                double nodeResizeRatio,
//                                                double cellResizeRatio ) :
//  Mesh( ndims, mesh_properties::mesh_of_cell_type[ CellType ], blockId,
// partId)
//// TODO: ???
////  m_node_coordinates( ndims, nodeCapacity, nodeResizeRatio ),
////  m_cell_connectivity( approxNumCellNodes( ndims, nodeCapacity,
////                                           cellNodeCapacity ),
////                       cellResizeRatio )
//{}

template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh( int ndims ) :
  Mesh( ndims, mesh_properties::mesh_of_cell_type[ CellType ], 0, 0 ),
  m_node_coordinates( new MeshCoordinates( ndims ) )
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh( int ndims,
                                                localIndex num_nodes,
                                                localIndex num_cells ) :
  Mesh( ndims, mesh_properties::mesh_of_cell_type[ CellType ], 0, 0 ),
  m_node_coordinates( new MeshCoordinates( ndims ) )
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::addCell( const localIndex* cell,
                                            int cell_type ) {
  SLIC_ASSERT( cell != AXOM_NULLPTR );
  SLIC_ASSERT( cell_type > MINT_UNDEFINED_CELL );
  SLIC_ASSERT( cell_type <= MINT_MIXED_CELL );
  SLIC_ASSERT( (CellType == MINT_MIXED_CELL) ? true : CellType == cell_type );

  m_cell_connectivity->addCell( cell, cell_type );

  this->setCellDataSize( getNumberOfCells() );
  // this->setFaceDataSize( WHAT???? )
  // this->setEdgeDataSize( WHAT???? )
}

//------------------------------------------------------------------------------
template < int CellType >
inline void UnstructuredMesh< CellType >::addNode( double x )
{
  SLIC_ASSERT( this->getDimension() == 1 );
  m_node_coordinates->addPoint( x );
  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline void UnstructuredMesh< CellType >::addNode( double x, double y )
{
  SLIC_ASSERT( this->getDimension() == 2 );
  m_node_coordinates->addPoint( x, y );
  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::addNode( double x, double y, double z ) {
  SLIC_ASSERT( this->getDimension() == 3 );
  m_node_coordinates->addPoint( x, y, z );
  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline void UnstructuredMesh< CellType >::addNode( const double* node ) {
  SLIC_ASSERT( node != AXOM_NULLPTR );

  if ( this->getDimension() == 1 )
  {
    m_node_coordinates->addPoint( node[0] );
  }
  else if ( this->getDimension() == 2 )
  {
    m_node_coordinates->addPoint( node[0], node[1] );
  }
  else
  {
    SLIC_ASSERT( this->getDimension() == 3 );
    m_node_coordinates->addPoint( node[0], node[1], node[2] );
  }

  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
const double* UnstructuredMesh< CellType >::getMeshCoordinateArray(int dim)
const {
  SLIC_ASSERT( dim < this->getDimension() );
  return m_node_coordinates->getCoordinateArray( dim );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
const localIndex* UnstructuredMesh< CellType >::getCell( localIndex cellIdx )
const
{
  SLIC_ASSERT(  cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
  return (*m_cell_connectivity)[ cellIdx ];
}

} /* namespace mint */
} /* namespace axom */

#endif /* UNSTRUCTUREDMESH_HXX_ */
