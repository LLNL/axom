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
#include "mint/MeshTypes.hpp"
#include "mint/CellTypes.hpp"
#include "mint/CellConnectivity.hpp"
#include "mint/MeshCoordinates.hpp"
#include "mint/Array.hpp"
#include "mint/config.hpp"

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
//  UnstructuredMesh( int ndims, IndexType nodeCapacity=100,
//                             IndexType cellNodeCapacity=0,
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
//                    IndexType nodeCapacity=100, IndexType
// cellNodeCapacity=0,
//                    double nodeResizeRatio=2.0, double cellResizeRatio=2.0 );

  UnstructuredMesh( int ndims );

  UnstructuredMesh( int ndims, IndexType num_nodes, IndexType num_cells );

  /*!
   * \brief Destructor.
   */
  virtual ~UnstructuredMesh() { }

  /*!
   * \brief Get the number of nodes that can be stored.
   * \return The capacity of the node array.
   */
  constexpr IndexType getNodeCapacity() const
  { return m_node_coordinates->capacity(); };

  constexpr double getNodeResizeRatio() const
  { return m_node_coordinates->getResizeRatio(); }


  void setNodeResizeRatio( double ratio )
  {
    m_node_coordinates->setResizeRatio( ratio );
    this->setNodeDataResizeRatio( ratio );
  }

  /*!
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return nnodes the number of nodes for the given cell.
   * \pre nnodes >= 1.
   */
  constexpr int getNumberOfCellNodes( IndexType cellIdx ) const
  { return m_cell_connectivity->getNumberOfNodes( cellIdx ); }

  /*!
   * \brief Get the number of nodes that can be stored in the cell
     connectivity.
   * \return the number of nodes that can be stored in the cell connectivity.
   */
  constexpr IndexType getCellNodeCapacity() const
  { return m_cell_connectivity->getCapacity(); }

  /*!
   * \brief Set the number of nodes that can be stored in the cell
     connectivity.
   * \param capacity the number of nodes that can be stored in the cell
   *  connectivity.
   */
  // void setCellNodeCapacity( IndexType capacity )
  // { m_cell_connectivity.setCellCapacity( capacity ); }


  constexpr double getCellResizeRatio() const
  { return m_cell_connectivity->getResizeRatio(); }


  void setCellResizeRatio( double ratio )
  {
    m_cell_connectivity->setResizeRatio( ratio );
    this->setCellDataResizeRatio( ratio );
    this->setFaceDataResizeRatio( ratio );
    this->setEdgeDataResizeRatio( ratio );
  }

  /*!
   * \brief Adds a new cell in the mesh
   * \param [in] cell pointer to the connectivity of the cell.
   * \param [in] cell_type the type of cell to Add.
   * \pre cell != AXOM_NULLPTR.
   */
  void addCell( const IndexType* cell, int cell_type );

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
   * \brief Returns pointer to the connectivity array of the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \return cell_ptr pointer to the connectivity array of the cell.
   * \pre cellIdx >= 0 && cellIdx < this->getNumberOfCells()
   * \post cell_ptr != AXOM_NULLPTR.
   */
  const IndexType* getCell( IndexType cellIdx ) const;
  IndexType* getCell( IndexType cellIdx );

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
  static IndexType approxNumCellNodes( int ndims, IndexType nodeCapacity,
                                        IndexType cellNodeCapacity )
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
    return std::pow( linear_num_cells, ndims )*cell_info[ CellType ].num_nodes;
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

// ------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh( int ndims ) :
  Mesh( ndims, mint::UNSTRUCTURED_MESH, 0, 0 ),
  m_node_coordinates( new MeshCoordinates( ndims ) ),
  m_cell_connectivity( new CellConnectivity< CellType >() )
{}

//------------------------------------------------------------------------------
// template < int CellType >
// UnstructuredMesh< CellType >::UnstructuredMesh( int ndims,
//                                                IndexType nodeCapacity,
//                                                IndexType cellNodeCapacity,
//                                                double nodeResizeRatio,
//                                                double cellResizeRatio ) :
//  Mesh( ndims, mesh_properties::mesh_of_cell_type[ CellType ], 0, 0 ),
//  m_node_coordinates( ndims, nodeCapacity, nodeResizeRatio ),
//  m_cell_connectivity( approxNumCellNodes( ndims, nodeCapacity,
//                                           cellNodeCapacity ), cellResizeRatio )
// {}

//------------------------------------------------------------------------------
// template < int CellType >
// UnstructuredMesh< CellType >::UnstructuredMesh( int ndims, int blockId,
//                                                int partId,
//                                                IndexType nodeCapacity,
//                                                IndexType cellNodeCapacity,
//                                                double nodeResizeRatio,
//                                                double cellResizeRatio ) :
//  Mesh( ndims, mesh_properties::mesh_of_cell_type[ CellType ], blockId,
// partId)
// TODO: ???
//  m_node_coordinates( ndims, nodeCapacity, nodeResizeRatio ),
//  m_cell_connectivity( approxNumCellNodes( ndims, nodeCapacity,
//                                           cellNodeCapacity ),
//                       cellResizeRatio )
// {}

// template < int CellType >
// UnstructuredMesh< CellType >::UnstructuredMesh( int ndims ) :
//   Mesh( ndims, mesh_properties::mesh_of_cell_type[ CellType ], 0, 0 ),
//   m_node_coordinates( new MeshCoordinates( ndims ) )
// {
//   TODO: implement this
// }

//------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh( int ndims,
                                                IndexType num_nodes,
                                                IndexType num_cells ) :
  Mesh( ndims, mint::UNSTRUCTURED_MESH, 0, 0 ),
  m_node_coordinates( new MeshCoordinates( num_nodes ) ),
  m_cell_connectivity( new CellConnectivity< CellType >( num_cells ) )
{
  // TODO: implement this
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::addCell( const IndexType* cell,
                                            int cell_type ) {
  SLIC_ASSERT( cell != AXOM_NULLPTR );
  SLIC_ASSERT( cell_type > mint::UNDEFINED_CELL );

  m_cell_connectivity->addCell( cell, cell_type );

// TODO: revisit this later
//  this->setCellDataSize( getNumberOfCells() );
  // this->setFaceDataSize( WHAT???? )
  // this->setEdgeDataSize( WHAT???? )
}

//------------------------------------------------------------------------------
template < int CellType >
inline void UnstructuredMesh< CellType >::addNode( double x )
{
  SLIC_ASSERT( this->getDimension() == 1 );
  m_node_coordinates->append( &x );
// TODO: revisit this later
//  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline void UnstructuredMesh< CellType >::addNode( double x, double y )
{
  SLIC_ASSERT( this->getDimension() == 2 );
  double xx[2] = { x, y };
  m_node_coordinates->append( xx );
// TODO: revisit this later
//  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::addNode( double x, double y, double z ) {
  SLIC_ASSERT( this->getDimension() == 3 );
  double xx[3] = {x, y, z};
  m_node_coordinates->append( xx );
// TODO: revisit this later
//  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline void UnstructuredMesh< CellType >::addNode( const double* node ) {
  SLIC_ASSERT( node != AXOM_NULLPTR );

  m_node_coordinates->append( node );
// TODO: revisiti this later
//  this->setNodeDataSize( getNumberOfNodes() );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
const IndexType* UnstructuredMesh< CellType >::getCell( IndexType cellIdx )
const
{
  // TODO: implement this
  return AXOM_NULLPTR;
//  SLIC_ASSERT(  cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
//  return (*m_cell_connectivity)[ cellIdx ];
}

//------------------------------------------------------------------------------
template < int CellType >
inline IndexType* UnstructuredMesh< CellType >::getCell( IndexType cellIdx )
{
  // TODO: implement this
  return AXOM_NULLPTR;
//  SLIC_ASSERT(  cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
//  return (*m_cell_connectivity)[ cellIdx ];
}

} /* namespace mint */
} /* namespace axom */

#endif /* UNSTRUCTUREDMESH_HXX_ */
