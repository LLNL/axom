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

#ifndef CELLCONNECTIVITY_HXX_
#define CELLCONNECTIVITY_HXX_

#include <iostream>

#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "mint/CellType.hpp"
#include "mint/Array.hpp"
#include "mint/config.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <cstring> // for memcpy()
#include <vector>  // for STL vector

namespace axom
{
namespace mint
{

template < int cell_type >
class CellConnectivity
{
public:

  /*!
   * \brief Default constructor.
   */
  CellConnectivity( IndexType capacity=100, double resize_ratio=2.0 ) :
    m_stride( cell::num_nodes[ cell_type ] ),
    m_connectivity( capacity, 0, cell::num_nodes[ cell_type ] )
  { m_connectivity.setResizeRatio( resize_ratio ); };

  /*!
   * \brief Destructor.
   */
  virtual ~CellConnectivity()
  {}

  /*!
   * \brief Checks if this CellConnecitivity instance has mixed cells.
   * \return status true if mixed cells else false.
   */
  bool hasMixedCellTypes() const
  { return false; }

  /*!
   * \brief Checks if this CellConnectivity instance is empty()
   * \return status true iff empty, else, false.
   */
  bool empty() const
  { return m_connectivity.size() == 0; }

  /*!
   * \brief Returns the total number of cells.
   * \return ncells number of cells.
   * \post ncells >= 0
   */
  IndexType getNumberOfCells() const
  { return m_connectivity.size(); }

  /*!
   * \brief Returns the number of nodes of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return nnodes number of nodes for the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post nnodes >= 0
   */
  IndexType getNumberOfNodes( IndexType AXOM_NOT_USED(cellIdx) ) const
  { return m_stride; };

  /*!
   * \brief Returns the cell type of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return ctype the cell type.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post ctype >= mint::VERTEX && ctype < mint::NUM_CELL_TYPES.
   */
  int getCellType( IndexType AXOM_NOT_USED(cellIdx) ) const
  { return cell_type; }

  /*!
   * \brief Access operator for the connectivity of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return cell_ptr pointer to the connectivity of the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post cell_ptr != AXOM_NULLPTR.
   */
  const IndexType* operator[]( IndexType cellIdx ) const {
    SLIC_ASSERT( ( cellIdx >= 0 ) && ( cellIdx < this->getNumberOfCells() ) );
    return m_connectivity.getData() + cellIdx * m_stride;
  }

  /*!
   * \brief Adds a cell in to this cell connectivity instance.
   * \param [in] cell array pointer to the connectivity of the cell to add.
   * \param [in] type the cell type.
   * \note type is only used for mixed cell connectivity.
   */
  void addCell( const IndexType* cell, int AXOM_NOT_USED(type) ) {
    SLIC_ASSERT( cell != AXOM_NULLPTR );
    m_connectivity.append( cell, m_stride );
  }

  /*!
   * \brief Sets the cell associated with the given cell index.
   * \param [in] cellIdx the index of the cell to set.
   * \param [in] cell array pointer to the source cell connectivity.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \pre cell != AXOM_NULLPTR
   */
  void setCell( IndexType cellIdx, const IndexType* cell ) {
    SLIC_ASSERT( ( cellIdx >= 0 ) && ( cellIdx < this->getNumberOfCells() ) );
    SLIC_ASSERT( cell != AXOM_NULLPTR );
    m_connectivity.set( cell, m_stride, cellIdx * m_stride );
  }

  /*!
   * \brief Get the maximum number of points that can currently be held.
   * \return N the capacity of m_coordinates.
   */
  IndexType getCapacity() const
  { return m_connectivity.getCapacity(); }


  void reserve( IndexType capacity )
  { m_connectivity.reserve( capacity ); }

  /*!
   * \brief Returns the number of points in this CellConnectivity instance.
   * \return npoint the number points in this CellConnectivity instance.
   */
  IndexType getTotalNumberOfNodes() const
  { return m_connectivity.size() * m_connectivity.getNumComponents(); }


  double getResizeRatio() const
  { return m_connectivity.getResizeRatio(); }


  void setResizeRatio( double ratio )
  { m_connectivity.setResizeRatio( ratio ); }

private:

  int m_stride;                          /*!< stride */
  Array< IndexType > m_connectivity;   /*!< connectivity array */

  CellConnectivity( const CellConnectivity& );
  CellConnectivity& operator=(const CellConnectivity& );
};

//------------------------------------------------------------------------------
//          Specialization of CellConnectivity for Mixed Elements
//------------------------------------------------------------------------------

template <>
class CellConnectivity< MINT_MIXED_CELL >
{
public:

  /*!
   * \brief Default constructor.
   */
  CellConnectivity( IndexType capacity=100, double resize_ratio=2.0 ) :
    m_num_cells(0),
    m_offset( capacity + 1, 0, 1 ),
    m_connectivity( capacity, 0, 1 ),
    m_cell_type( capacity, 0, 1 )
  {
    //
    // TO DO: Need capacity for cells and nodes.
    //
    m_offset.setResizeRatio( resize_ratio );
    m_connectivity.setResizeRatio( resize_ratio );
    m_cell_type.setResizeRatio( resize_ratio );
  };

  /*!
   * \brief Destructor.
   */
  virtual ~CellConnectivity()
  {};

  /*!
   * \brief Checks if this CellConnecitivity instance has mixed cells.
   * \return status true if mixed cells else false.
   */
  bool hasMixedCellTypes() const
  { return true; }

  /*!
   * \brief Checks if this CellConnectivity instance is empty()
   * \return status true iff empty, else, false.
   */
  bool empty() const
  { return this->getNumberOfCells() == 0; };

  /*!
   * \brief Returns the total number of cells.
   * \return ncells number of cells.
   * \post ncells >= 0
   */
  IndexType getNumberOfCells() const
  { return m_num_cells; };

  /*!
   * \brief Returns the number of nodes of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return nnodes number of nodes for the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post nnodes >= 0
   */
  IndexType getNumberOfNodes( IndexType cellIdx ) const {
    SLIC_ASSERT( cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
    return m_offset( cellIdx + 1 ) - m_offset( cellIdx );
  };

  /*!
   * \brief Returns the cell type of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return ctype the cell type.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post ctype >= mint::VERTEX && ctype < mint::NUM_CELL_TYPES.
   */
  int getCellType( IndexType cellIdx ) const {
    SLIC_ASSERT( cellIdx >= 0 && cellIdx < this->getNumberOfCells() );
    return m_cell_type( cellIdx );
  }

  /*!
   * \brief Access operator for the connectivity of the given cell.
   * \param [in] cellIdx the index of the cell in question.
   * \return cell_ptr pointer to the connectivity of the given cell.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \post cell_ptr != AXOM_NULLPTR.
   */
  const IndexType* operator[]( IndexType cellIdx ) const {
    SLIC_ASSERT( ( cellIdx >= 0 ) && ( cellIdx < this->getNumberOfCells() ) );
    return m_connectivity.getData() + m_offset( cellIdx );
  };

  /*!
   * \brief Addss a cell in to this cell connectivity instance.
   * \param [in] cell array pointer to the connectivity of the cell to add.
   * \param [in] type the cell type.
   * \param [in] nnodes the number of nodes in the cell.
   * \note type is only used for mixed cell connectivity.
   * \pre cell != AXOM_NULLPTR .
   */
  void addCell( const IndexType* cell, int type ) {
    int num_nodes = cell::num_nodes[ type ];

    /* STEP 0: get the last cell index before adding the new cell. */
    IndexType new_cell_id  = this->getNumberOfCells();

    if ( this->empty() )
    {
      m_offset.append( 0 );
    }

    /* STEP 2: update the offsets array. */
    const IndexType offset = m_offset( new_cell_id );
    m_offset.append( offset + num_nodes );

    /* STEP 3: update the cell connectivity. */
    m_connectivity.append( cell, num_nodes );

    /* STEP 4: update the cell types. */
    m_cell_type.append( type );

    m_num_cells++;
  }

  /*!
   * \brief Sets the cell associated with the given cell index.
   * \param [in] cellIdx the index of the cell to set.
   * \param [in] cell array pointer to the source cell connectivity.
   * \pre cellIdx >= 0 && cellIdx < ncells
   * \pre cell != AXOM_NULLPTR
   */
  void setCell( IndexType cellIdx, const IndexType* cell )
  {
    SLIC_ASSERT( ( cellIdx >= 0 ) && ( cellIdx < this->getNumberOfCells() ) );
    SLIC_ASSERT( cell != AXOM_NULLPTR );

    /* STEP 0: get the number of nodes for the given cell type. */
    const IndexType nnodes = this->getNumberOfNodes( cellIdx );

    /* STEP 1: get to/from pointers. */
    const IndexType offset = m_offset( cellIdx );

    m_connectivity.set( cell, nnodes, offset );
  }

  /*!
   * \brief Get the maximum number of points that can currently be held.
   * \return N the capacity of m_coordinates.
   */
  IndexType getCapacity() const
  { return m_connectivity.getCapacity(); }


  void reserve( IndexType capacity )
  {
    //
    // TO DO: Need capacity for cells and nodes.
    //
    m_offset.reserve( capacity + 1 );
    m_connectivity.reserve( capacity );
    m_cell_type.reserve( capacity );
  }

  /*!
   * \brief Returns the number of points in this CellConnectivity instance.
   * \return npoint the number points in this CellConnectivity instance.
   */
  IndexType getSize() const
  { return m_connectivity.size(); }


  double getResizeRatio() const
  { return m_connectivity.getResizeRatio(); }


  void setResizeRatio( double ratio )
  {
    m_offset.setResizeRatio( ratio );
    m_connectivity.setResizeRatio( ratio );
    m_cell_type.setResizeRatio( ratio );
  }

private:
  IndexType m_num_cells;
  Array< IndexType > m_offset;
  Array< IndexType > m_connectivity;
  Array< unsigned char > m_cell_type;

  CellConnectivity( const CellConnectivity& );
  CellConnectivity& operator=( const CellConnectivity& );
};

} /* namespace mint */
} /* namespace axom */

#endif /* CELLCONNECTIVITY_HXX_ */
