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

#ifndef MINT_CELLTYPES_HPP_
#define MINT_CELLTYPES_HPP_

#include "mint/config.hpp"

namespace axom
{
namespace mint
{

constexpr int MAX_NUM_NODES = 27;

enum class Topology : unsigned char
{
  SINGLE,
  MIXED
};


enum class ConnectivityType : unsigned char
{
  NO_INDIRECTION,
  INDIRECTION,
  TYPED_INDIRECTION
};


template < Topology TOPO >
struct topology_traits {};

template <>
struct topology_traits< Topology::SINGLE > 
{
  constexpr static ConnectivityType cell_connec = ConnectivityType::NO_INDIRECTION;
};

template <>
struct topology_traits< Topology::MIXED > 
{
  constexpr static ConnectivityType cell_connec = ConnectivityType::TYPED_INDIRECTION;
};


/*!
 * \brief Enumerates all cell types supported by Mint
 */
using CellType = signed char;
enum CellTypes : CellType
{
  UNDEFINED_CELL = -1,    ///< UNDEFINED

  VERTEX,                 ///< VERTEX
  SEGMENT,                ///< LINE_SEGMENT

  TRIANGLE,               ///< LINEAR_TRIANGLE
  QUAD,                   ///< LINEAR_QUAD
  TET,                    ///< LINEAR_TET
  HEX,                    ///< LINEAR_HEX
  PRISM,                  ///< LINEAR_PRISM
  PYRAMID,                ///< LINEAR_PYRAMID

  QUAD9,                  ///< QUADRATIC QUAD
  HEX27,                  ///< QUADRATIC HEX

  NUM_CELL_TYPES          ///<  total number of cell types
};

/*!
 * \def REGISTER_CELL_INFO( MINT_CELL_TYPE, MINT_NAME, BP_NAME, VTK_TYPE, N )
 *
 * \brief Convenience macro used to register information about a cell type.
 *
 * \param MINT_CELL_TYPE the mint cell type, e.g., mint::QUAD, mint::HEX, etc.
 * \param MINT_NAME the associated mint name for the cell type.
 * \param BP_NAME the associated name in the mesh blueprint.
 * \param VTK_TYPE the corresponding VTK type.
 * \param N the number of nodes that the cell has.
 */
#define REGISTER_CELL_INFO( MINT_CELL_TYPE, MINT_NAME, BP_NAME, VTK_TYPE, N ) \
  namespace internal                                                          \
  {                                                                           \
    static constexpr CellInfo MINT_CELL_TYPE##_INFO =                         \
    {                                                                         \
      MINT_CELL_TYPE,                                                         \
      MINT_NAME,                                                              \
      BP_NAME,                                                                \
      VTK_TYPE,                                                               \
      N                                                                       \
    };                                                                        \
  }

/*!
 * \def CELL_INFO( MINT_CELL_TYPE )
 *
 * \brief Convenience macro used to access a registered CellInfo struct for
 *  the specified mint cell type.
 *
 * \param MINT_CELL_TYPE the mint cell type, e.g., mint::QUAD, mint::HEX, etc.
 */
#define CELL_INFO( MINT_CELL_TYPE ) internal::MINT_CELL_TYPE##_INFO

/*!
 * \struct CellInfo
 *
 * \brief Holds information associated with a given cell type.
 */
typedef struct
{
   CellType cell_type;          /*!< cell type, .e.g, mint::QUAD, mint::HEX */
   const char* name;            /*!< the name associated with the cell */
   const char* blueprint_name;  /*!< corresponding mesh blueprint name */
   int vtk_type;                /*!< corresponding vtk_type */
   int num_nodes;               /*!< number of nodes for the given cell */
} CellInfo;

// Cell Info registration
REGISTER_CELL_INFO( VERTEX, "VERTEX", "point", 1, 1 );
REGISTER_CELL_INFO( SEGMENT, "SEGMENT", "line", 3, 2 );

REGISTER_CELL_INFO( TRIANGLE, "TRIANGLE", "tri", 5, 3 );
REGISTER_CELL_INFO( QUAD, "QUAD", "quad", 9, 4 );
REGISTER_CELL_INFO( TET, "TET", "tet", 10, 4 );
REGISTER_CELL_INFO( HEX, "HEX", "hex", 12, 8 );
REGISTER_CELL_INFO( PRISM, "PRISM", "prism-no-bp", 13, 6 );
REGISTER_CELL_INFO( PYRAMID, "PYRAMID", "pyramid-no-bp", 14, 5 );

REGISTER_CELL_INFO( QUAD9, "QUAD9", "quad9-no-bp", 28, 9 );
REGISTER_CELL_INFO( HEX27, "HEX27", "hex27-no-bp", 29, 27 );

/*!
 * \brief Array of CellInfo corresponding to each cell type
 * \note The order at which CellInfo for each type is added has to match
 *  the order of the cell types in the CellTypes enum above.
 */
static constexpr CellInfo cell_info[ NUM_CELL_TYPES ] = {
  CELL_INFO( VERTEX ),
  CELL_INFO( SEGMENT ),
  CELL_INFO( TRIANGLE ),
  CELL_INFO( QUAD ),
  CELL_INFO( TET ),
  CELL_INFO( HEX ),
  CELL_INFO( PRISM ),
  CELL_INFO( PYRAMID ),
  CELL_INFO( QUAD9 ),
  CELL_INFO( HEX27 )
};


} /* namespace mint */
} /* namespace axom */

#endif /* MINT_CellTypes_HXX_ */
