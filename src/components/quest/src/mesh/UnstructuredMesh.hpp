/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file UnstructuredMesh.hxx
 *
 * \date Sep 13, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef UNSTRUCTUREDMESH_HXX_
#define UNSTRUCTUREDMESH_HXX_

#include "quest/Mesh.hpp"
#include "quest/MeshType.hpp"
#include "quest/CellType.hpp"
#include "quest/CellConnectivity.hpp"
#include "quest/MeshCoordinates.hpp"
#include "quest/Field.hpp"
#include "quest/FieldTypes.hpp"
#include "quest/FieldData.hpp"

#include "slic/slic.hpp"

#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"

#include <cstring> // for memcpy()
#include <fstream> // for fstream

namespace meshtk {

template < int CellType >
class UnstructuredMesh : public Mesh
{
public:

  /*!
   *****************************************************************************
   * \brief Custom constructor. Creates an unstructured mesh of given dimension.
   * \param [in] ndims the number of dimension.
   *****************************************************************************
   */
  explicit UnstructuredMesh(int ndims);

  /*!
   *****************************************************************************
   * \brief Creates an unstructured mesh of given dimension that is identified
   *  by the supplied blockId and partitionId pair.
   * \param [in] ndims the number of dimension.
   * \param [in] blockId the block ID of the mesh
   * \param [in] partId the partition ID of the mesh
   *****************************************************************************
   */
  UnstructuredMesh(int ndims, int blockId, int partId);

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  virtual ~UnstructuredMesh();

  /// \name Virtual API
  /// @{

  /*!
   *****************************************************************************
   * \brief Returns the total number of nodes in the mesh.
   * \return numNodes the total number of nodes.
   * \post numNodes >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   *****************************************************************************
   */
  virtual int getMeshNumberOfNodes() const
    { return this->getNumberOfNodes(); };

  /*!
   *****************************************************************************
   * \brief Returns the total number of cells in the mesh.
   * \return numCells the total number of cells.
   * \post numCells >= 0
   * \warning This is a virtual method -- do not call inside a loop.
   *****************************************************************************
   */
  virtual int getMeshNumberOfCells() const
    { return this->getNumberOfCells(); };

  /*!
   *****************************************************************************
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return numCellNodes the number of nodes in the given cell.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   *****************************************************************************
   */
  virtual int getMeshNumberOfCellNodes( int cellIdx ) const
    { return this->getNumberOfCellNodes( cellIdx ); };

  /*!
   *****************************************************************************
   * \brief Returns the cell connectivity of the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \param [out] cell user-supplied buffer to store cell connectivity info.
   * \note cell must have sufficient size to hold the connectivity information.
   * \pre cellIdx >= 0 && cellIdx < this->getMeshNumberOfCells()
   * \pre cell != ATK_NULLPTR.
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   *****************************************************************************
   */
  virtual void getMeshCell( int cellIdx, int* cell ) const
  {
    const int numNodes = this->getNumberOfCellNodes( cellIdx );
    const int* myCell  = this->getCell( cellIdx );
    memcpy( cell, myCell, numNodes*sizeof(int) );
  }

  /*!
   *****************************************************************************
   * \brief Returns the cell type of the cell associated with the given Id.
   * \param [in] cellIdx the index of the cell in query.
   * \return cellType the cell type of the cell at the given index.
   *****************************************************************************
   */
  virtual int getMeshCellType( int cellIdx ) const
    { return m_cell_connectivity->getCellType( cellIdx ); }

  /*!
   *****************************************************************************
   * \brief Returns the coordinates of the given node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] coordinates user-supplied buffer to store the node coordinates.
   * \pre nodeIdx >= && nodeIdx < this->getMeshNumberOfNodes()
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   *****************************************************************************
   */
  virtual void getMeshNode( int nodeIdx, double* coordinates ) const
  {
    for ( int i=0; i < m_ndims; ++i ) {
      const double* coord = this->getMeshCoordinateArray( i );
      coordinates[ i ] = coord[ nodeIdx ];
    } // END for all dimensions
  }

  /*!
   *****************************************************************************
   * \brief Returns the coordinate of a mesh node.
   * \param [in] nodeIdx the index of the node in query.
   * \param [in] dim the dimension of the coordinate to return, e.g., x, y or z
   * \return c the coordinate value of the node at
   * \pre dim >= 0 && dim < m_ndims
   * \warning this is a virtual method, downcast to the derived class and use
   *  the non-virtual API instead to avoid the overhead of a virtual call.
   *****************************************************************************
   */
  virtual double getMeshNodeCoordinate( int nodeIdx, int dim ) const
  {
    SLIC_ASSERT( dim >= 0 && dim < m_ndims );
    return this->getMeshCoordinateArray( dim )[ nodeIdx ];
  }


  /// @}

  /*!
   *****************************************************************************
   * \brief Returns the total number of nodes in the mesh.
   * \return N the total number of nodes in the mesh.
   * \pre m_node_coordinates != ATK_NULLPTR.
   * \post N >= 0.
   *****************************************************************************
   */
  int getNumberOfNodes() const
    { return m_node_coordinates->getNumberOfPoints(); };

  /*!
   *****************************************************************************
   * \brief Returns the total number cells in the mesh.
   * \return N the total number of cells in the mesh.
   * \pre m_cell_connectivity != ATK_NULLPTR.
   * \post N >= 0.
   *****************************************************************************
   */
  int getNumberOfCells() const
    { return m_cell_connectivity->getNumberOfCells(); };

  /*!
   *****************************************************************************
   * \brief Returns the number of nodes for the given cell.
   * \param cellIdx the index of the cell in query.
   * \return nnodes the number of nodes for the given cell.
   * \pre nnodes >= 1.
   *****************************************************************************
   */
  int getNumberOfCellNodes( int cellIdx ) const
    { return m_cell_connectivity->getNumberOfNodes( cellIdx ); };

  /*!
   *****************************************************************************
   * \brief Inserts a new cell in the mesh
   * \param [in] cell pointer to the connectivity of the cell.
   * \param [in] cell_type the type of cell to insert.
   * \param [in] numNodes number of nodes for the cell.
   * \pre cell != ATK_NULLPTR.
   *****************************************************************************
   */
  void insertCell( const int* cell, int cell_type, int numNodes );

  /*!
   *****************************************************************************
   * \brief Inserts a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \param [in] y the y--coordinate of the new mesh node.
   * \pre m_ndims == 2.
   *****************************************************************************
   */
  void insertNode( double x, double y );

  /*!
   *****************************************************************************
   * \brief Inserts a new node in the mesh.
   * \param [in] x the x--coordinate of the new mesh node.
   * \param [in] y the y--coordinate of the new mesh node.
   * \param [in] z the z--coordinate of the new mesh node.
   * \pre m_ndims == 3.
   *****************************************************************************
   */
  void insertNode( double x, double y, double z );

  /*!
   *****************************************************************************
   * \brief Returns pointer to coordinates array for the requested dimension.
   * \param [in] dim the requested dimension.
   * \return ptr pointer to the coordinates array for the given dimension.
   * \pre m_coordinates != ATK_NULLPTR.
   * \pre dim < m_ndims
   * \post ptr != ATK_NULLPTR.
   *****************************************************************************
   */
  const double* getMeshCoordinateArray( int dim ) const;

  /*!
   *****************************************************************************
   * \brief Returns pointer to the connectivity array of the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \return cell_ptr pointer to the connectivity array of the cell.
   * \pre m_cell_connectivity != ATK_NULLPTR
   * \pre cellIdx >= 0 && cellIdx < this->getNumberOfCells()
   * \post cell_ptr != ATK_NULLPTR.
   *****************************************************************************
   */
  const int* getCell( int cellIdx ) const;

  /*!
   *****************************************************************************
   * \brief Writes a VTK file corresponding to this mesh instance in the VTK
   *  Legacy ASCII format that can be visualized with VisIt or ParaView.
   * \param [in] vtkFileName the name of the VTK file.
   * \note This method is primarily intended for debugging.
   *****************************************************************************
   */
  void toVtkFile(const std::string& vtkFileName);

private:

  /*!
   *****************************************************************************
   * \brief Default constructor.
   * \note Made private to prevent users from calling it.
   *****************************************************************************
   */
  UnstructuredMesh();


  MeshCoordinates*                      m_node_coordinates;
  CellConnectivity< int, CellType >*    m_cell_connectivity;

  DISABLE_COPY_AND_ASSIGNMENT(UnstructuredMesh);
};

} /* namespace meshtk */

//------------------------------------------------------------------------------
//      UnstructuredMesh Implementation
//------------------------------------------------------------------------------
namespace meshtk {

template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh() :
    Mesh(-1,-1,-1,-1),
    m_node_coordinates( ATK_NULLPTR ),
    m_cell_connectivity( ATK_NULLPTR )
{

}

//------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh(int ndims ) :
    Mesh(ndims, mesh_properties::mesh_of_cell_type[ CellType ], 0, 0 ),
    m_node_coordinates( new MeshCoordinates( ndims ) ),
    m_cell_connectivity( new CellConnectivity< int, CellType >() )
{

}

//------------------------------------------------------------------------------
template < int CellType >
UnstructuredMesh< CellType >::UnstructuredMesh(
    int ndims, int blockId, int partId ) :
    Mesh(ndims, mesh_properties::mesh_of_cell_type[CellType],blockId,partId ),
    m_node_coordinates( new MeshCoordinates( ndims ) ),
    m_cell_connectivity( new CellConnectivity< int, CellType >() )
{

}

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
void UnstructuredMesh< CellType >::insertCell(
    const int* cell, int cell_type, int num_nodes )
{
  SLIC_ASSERT( cell != ATK_NULLPTR );
  SLIC_ASSERT( m_cell_connectivity != ATK_NULLPTR );
  SLIC_ASSERT( cell::num_nodes[ cell_type ] == num_nodes );

  m_cell_connectivity->insertCell( cell, cell_type, num_nodes );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::insertNode( double x, double y )
{
  SLIC_ASSERT( m_ndims==2 );
  SLIC_ASSERT( m_node_coordinates != ATK_NULLPTR );
  m_node_coordinates->insertPoint( x, y );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::insertNode( double x, double y, double z )
{
  SLIC_ASSERT( m_ndims==3 );
  SLIC_ASSERT( m_node_coordinates != ATK_NULLPTR );
  m_node_coordinates->insertPoint( x, y, z );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
const double* UnstructuredMesh< CellType >::getMeshCoordinateArray(int dim) const
{
  SLIC_ASSERT( m_node_coordinates != ATK_NULLPTR );
  SLIC_ASSERT( dim < m_ndims );
  return m_node_coordinates->getCoordinateArray( dim );
}

//------------------------------------------------------------------------------
template < int CellType >
inline
const int* UnstructuredMesh< CellType >::getCell( int cellIdx ) const
{
  SLIC_ASSERT( m_cell_connectivity != ATK_NULLPTR );
  SLIC_ASSERT( cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
  return (*m_cell_connectivity)[ cellIdx ];
}

//------------------------------------------------------------------------------
template < int CellType >
inline
void UnstructuredMesh< CellType >::toVtkFile( const std::string& file )
{
  std::ofstream ofs;
  ofs.open( file.c_str() );

  // STEP 0: Write VTK header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " Unstructured Mesh\n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // STEP 1: Write mesh nodes
  ofs << "POINTS " << this->getNumberOfNodes() << " double\n";
  for ( int node=0; node < this->getNumberOfNodes(); ++node ) {

    ofs << m_node_coordinates->getCoordinate( node, 0 ) << " ";
    ofs << m_node_coordinates->getCoordinate( node, 1 ) << " ";
    if ( m_ndims==3 ) {

        ofs << m_node_coordinates->getCoordinate( node, 2 );

    } else {

        ofs << "0.0";
    }

    ofs << std::endl;

  } // END for all nodes

  // STEP 2: Write mesh cell connectivity
  const int ncells   = m_cell_connectivity->getNumberOfCells();

  int numCellNodes = 0;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ){
      numCellNodes += m_cell_connectivity->getNumberOfNodes( cellIdx ) + 1;
  }
  ofs << "CELLS " << ncells << " " << numCellNodes << std::endl;

  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {

    const int nnodes = m_cell_connectivity->getNumberOfNodes( cellIdx );
    const int* cell  = (*m_cell_connectivity)[ cellIdx ];

    ofs << nnodes << " ";
    for ( int i=0; i < nnodes; ++i ) {
      ofs << cell[ i ] << " ";
    } // END for all nodes
    ofs << std::endl;

  } // END for all cells

  // STEP 3: Write cell types
  ofs << "CELL_TYPES " << ncells << std::endl;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {
    int ctype    = m_cell_connectivity->getCellType( cellIdx );
    int vtk_type = cell::vtk_types[ ctype ];
    ofs << vtk_type << std::endl;
  } // END for all cells

  // STEP 4: Write Cell Data
  ofs << "CELL_DATA " << ncells << std::endl;
  FieldData* CD = this->getCellFieldData();
  for ( int f=0; f < CD->getNumberOfFields(); ++f ) {

      Field* field = CD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i ) {
             ofs << dataPtr[ i ] << std::endl;
          }
      }
  }

  // STEP 5: Write Point Data
  const int nnodes = this->getMeshNumberOfNodes();
  ofs << "POINT_DATA " << nnodes << std::endl;
  FieldData* PD = this->getNodeFieldData();
  for ( int f=0; f < PD->getNumberOfFields(); ++f ) {

      Field* field = PD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }
      }
  }

  ofs.close();
}

} /* namespace meshtk */

#endif /* UNSTRUCTUREDMESH_HXX_ */
