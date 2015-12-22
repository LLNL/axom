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
 * \file CellConnectivity.hxx
 *
 * \date Sep 12, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */


#ifndef CELLCONNECTIVITY_HXX_
#define CELLCONNECTIVITY_HXX_

#define NotUsed(x)

#include <iostream>

#include "quest/CellType.hpp"

#include "common/CommonTypes.hpp"

// C/C++ includes
#include <cassert> // for assert()
#include <cstring> // for memcpy()
#include <vector>  // for STL vector

namespace meshtk {


template < typename index_type, int cell_type >
class CellConnectivity
{
public:

  /*!
   *****************************************************************************
   * \brief Default constructor.
   *****************************************************************************
   */
  CellConnectivity() : m_stride( cell::num_nodes[ cell_type ] )
  {
    m_connectivity.reserve( 100*m_stride );
  };

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  virtual ~CellConnectivity() { m_connectivity.clear(); };

  /*!
   *****************************************************************************
   * \brief Checks if this CellConnecitivity instance has mixed cells.
   * \return status true if mixed cells else false.
   *****************************************************************************
   */
  bool hasMixedCellTypes() const { return false; }

  /*!
   *****************************************************************************
   * \brief Checks if this CellConnectivity instance is empty()
   * \return status true iff empty, else, false.
   *****************************************************************************
   */
  bool empty() const { return (this->getNumberOfCells()==0); };

  /*!
   *****************************************************************************
   * \brief Returns the total number of cells.
   * \return ncells number of cells.
   * \post ncells >= 0
   *****************************************************************************
   */
  int getNumberOfCells() const { return( m_connectivity.size()/m_stride ); };

  /*!
   *****************************************************************************
   * \brief Returns the number of nodes of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return nnodes number of nodes for the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post nnodes >= 0
   *****************************************************************************
   */
  int getNumberOfNodes( int NotUsed(cellIdx) ) const { return m_stride; };

  /*!
   *****************************************************************************
   * \brief Returns the cell type of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return ctype the cell type.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post ctype >= meshtk::VERTEX && ctype < meshtk::NUM_CELL_TYPES.
   *****************************************************************************
   */
  int getCellType( int NotUsed(cellIdx) ) const { return cell_type; }

  /*!
   *****************************************************************************
   * \brief Access operator for the connectivity of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return cell_ptr pointer to the connectivity of the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post cell_ptr != ATK_NULLPTR.
   *****************************************************************************
   */
  const index_type* operator[]( int cellIdx )
  {
    assert( (cellIdx >= 0) && (cellIdx < this->getNumberOfCells() ) );
    return &m_connectivity[ cellIdx*m_stride ];
  };

  /*!
   *****************************************************************************
   * \brief Inserts a cell in to this cell connectivity instance.
   * \param [in] cell array pointer to the connectivity of the cell to insert.
   * \param [in] type the cell type.
   * \note type is only used for mixed cell connectivity.
   *****************************************************************************
   */
  void insertCell( const index_type* cell,
                   int NotUsed(type),
                   int NotUsed(nnodes) )
  {
    assert( cell != ATK_NULLPTR );

    for ( int i=0; i < m_stride; ++i ) {
      m_connectivity.push_back( cell[ i ] );
    }

  }

  /*!
   *****************************************************************************
   * \brief Sets the cell associated with the given cell index.
   * \param [in] cellIdx the index of the cell to set.
   * \param [in] cell array pointer to the source cell connectivity.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \pre cell != ATK_NULLPTR
   *****************************************************************************
   */
  void setCell( int cellIdx, const index_type* cell )
  {
    assert( (cellIdx >= 0) && (cellIdx < this->getNumberOfCells()) );
    assert( cell != ATK_NULLPTR );

    index_type* to         = &m_connectivity[ cellIdx*m_stride ];
    const index_type* from = cell;
    memcpy( to, from, m_stride*sizeof( index_type )  );
  }

private:
  int m_stride;                             /*!< stride */
  std::vector< index_type > m_connectivity; /*!< connectivity array */

  CellConnectivity( const CellConnectivity& );
  CellConnectivity& operator=(const CellConnectivity& );
};

//------------------------------------------------------------------------------
//          Specialization of CellConnectivity for Mixed Elements
//------------------------------------------------------------------------------

template < typename index_type >
class CellConnectivity< index_type, MIXED >
{
public:

  /*!
   *****************************************************************************
   * \brief Default constructor.
   *****************************************************************************
   */
  CellConnectivity() : m_num_cells(0)
  {
    m_offset.reserve( 101 );
    m_cell_type.reserve( 100 );
    m_connectivity.reserve( 100*cell::num_nodes[ MIXED ] );
  };

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  virtual ~CellConnectivity()
  {
    m_offset.clear();
    m_connectivity.clear();
  };

  /*!
   *****************************************************************************
   * \brief Checks if this CellConnecitivity instance has mixed cells.
   * \return status true if mixed cells else false.
   *****************************************************************************
   */
  bool hasMixedCellTypes() const { return true; }

  /*!
   *****************************************************************************
   * \brief Checks if this CellConnectivity instance is empty()
   * \return status true iff empty, else, false.
   *****************************************************************************
   */
  bool empty() const { return (this->getNumberOfCells()==0); };

  /*!
   *****************************************************************************
   * \brief Returns the total number of cells.
   * \return ncells number of cells.
   * \post ncells >= 0
   *****************************************************************************
   */
  int getNumberOfCells() const { return m_num_cells; };

  /*!
   *****************************************************************************
   * \brief Returns the number of nodes of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return nnodes number of nodes for the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post nnodes >= 0
   *****************************************************************************
   */
  int getNumberOfNodes( int cellIdx ) const
  {
    assert( cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
    return m_offset[ cellIdx+1 ] - m_offset[ cellIdx ];
  };

  /*!
   *****************************************************************************
   * \brief Returns the cell type of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return ctype the cell type.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post ctype >= meshtk::VERTEX && ctype < meshtk::NUM_CELL_TYPES.
   *****************************************************************************
   */
  int getCellType( int cellIdx ) const
  {
    assert( cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
    return m_cell_type[ cellIdx ];
  }

  /*!
   *****************************************************************************
   * \brief Access operator for the connectivity of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return cell_ptr pointer to the connectivity of the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post cell_ptr != ATK_NULLPTR.
   *****************************************************************************
   */
  const index_type* operator[]( int cellIdx )
  {
    assert( (cellIdx >= 0) && (cellIdx < this->getNumberOfCells() ) );
    return &m_connectivity[ m_offset[cellIdx] ];
  };

  /*!
   *****************************************************************************
   * \brief Inserts a cell in to this cell connectivity instance.
   * \param [in] cell array pointer to the connectivity of the cell to insert.
   * \param [in] type the cell type.
   * \param [in] nnodes the number of nodes in the cell.
   * \note type is only used for mixed cell connectivity.
   * \pre cell != ATK_NULLPTR .
   *****************************************************************************
   */
  void insertCell( const index_type* cell, int type, int nnodes )
  {
    assert( cell != ATK_NULLPTR );

    // STEP 0: get the last cell index before inserting the new cell
    int last_cell_id = this->getNumberOfCells()-1;
    int new_cell_id  = last_cell_id+1;

    if ( this->empty() ) {
      m_offset.push_back(0);
    }

    // STEP 2: update the offsets array
    const int offset = m_offset[ new_cell_id ];
    m_offset.push_back( offset+nnodes );


    // STEP 3: update the cell connectivity
    for ( int i=0; i < nnodes; ++i ) {
      m_connectivity.push_back( cell[ i ] );
    } // END for all nodes

    m_cell_type.push_back( type );

    ++m_num_cells;
  }

  /*!
   *****************************************************************************
   * \brief Sets the cell associated with the given cell index.
   * \param [in] cellIdx the index of the cell to set.
   * \param [in] cell array pointer to the source cell connectivity.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \pre cell != ATK_NULLPTR
   *****************************************************************************
   */
  void setCell( int cellIdx, const index_type* cell )
  {
    assert( (cellIdx >= 0) && (cellIdx < this->getNumberOfCells()) );
    assert( cell != ATK_NULLPTR );

    // STEP 0: get the number of nodes for the given cell type
    const int nnodes = this->getNumberOfNodes( cellIdx );

    // STEP 1: get to/from pointers
    const int offset       = m_offset[ cellIdx ];
    index_type* to         = &m_connectivity[ offset ];
    const index_type* from = cell;

    // STEP 2: bytecopy the connectivity information
    memcpy( to, from, nnodes*sizeof(index_type) );
  }

private:
  int m_num_cells;
  std::vector< int >        m_offset;
  std::vector< index_type > m_connectivity;
  std::vector < int >       m_cell_type;

  CellConnectivity( const CellConnectivity& );
  CellConnectivity& operator=(const CellConnectivity& );
};

/// \name CellConnectivity Type Definitions
/// \brief Defined here for convenience.
/// @{

typedef CellConnectivity< int, LINE >     LineCellConnectivity;
typedef CellConnectivity< int, LINEAR_QUAD >     QuadCellConnectivity;
typedef CellConnectivity< int, LINEAR_TRIANGLE > TriangleCellConnectivity;
typedef CellConnectivity< int, LINEAR_TET >      TetCellConnectivity;
typedef CellConnectivity< int, LINEAR_HEX >      HexCellConnectivity;

typedef CellConnectivity< int, MIXED  > MixedCellConnectivity;
/// @}

} /* namespace meshtk */

#endif /* CELLCONNECTIVITY_HXX_ */
