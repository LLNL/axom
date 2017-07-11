/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef CELLTYPE_HPP_
#define CELLTYPE_HPP_

#include <string>

/*!
 * \file
 *
 * \enum CellType
 * \brief Defines the basic cell types supported by mint.
 */
enum {
  MINT_UNDEFINED_CELL = -1, //!< UNDEFINED

  MINT_VERTEX,         ///< VERTEX
  MINT_SEGMENT,        ///< LINE_SEGMENT

  MINT_TRIANGLE,       ///< LINEAR_TRIANGLE
  MINT_QUAD,           ///< LINEAR_QUAD
  MINT_TET,            ///< LINEAR_TET
  MINT_HEX,            ///< LINEAR_HEX
  MINT_PRISM,          ///< LINEAR_PRISM
  MINT_PYRAMID,        ///< LINEAR_PYRAMID

  MINT_MIXED_CELL,     ///< MIXED
  MINT_NUM_CELL_TYPES  ///< NUM_CELL_TYPES
};

#define MINT_MAX_NUM_NODES 8

namespace axom {
namespace mint {

namespace cell {

static const int vtk_types[] = {
  1,    // VERTEX          -> VTK_VERTEX
  3,    // LINE_SEGMENT    -> VTK_LINE
  5,    // LINEAR_TRIANGLE -> VTK_TRIANGLE
  9,    // LINEAR_QUAD     -> VTK_QUAD
  10,   // LINEAR_TET      -> VTK_TET
  12,   // LINEAR_HEX      -> VTK_HEXAHEDRON
  13,   // LINEAR_PRISM    -> VTK_WEDGE
  14,   // LINEAR_PYRAMID  -> VTK_PYRAMID
  0     // MIXED           -> VTK_EMPTY_CELL
};

static const int num_nodes[] = {
  1,    // VERTEX
  2,    // LINE_SEGMENT
  3,    // LINEAR_TRIANGLE
  4,    // LINEAR_QUAD
  4,    // LINEAR_TET
  8,    // LINEAR_HEX
  6,    // LINEAR_PRISM
  5,    // LINEAR_PYRAMID

  5     // MIXED
};

static const std::string name[] = {
  "MINT_VERTEX",         //!< VERTEX
  "MINT_SEGMENT",        //!< LINE_SEGMENT
  "MINT_TRIANGLE",       //!< LINEAR_TRIANGLE
  "MINT_QUAD",           //!< LINEAR_QUAD
  "MINT_TET",            //!< LINEAR_TET
  "MINT_HEX",            //!< LINEAR_HEX
  "MINT_PRISM",          //!< LINEAR_PRISM
  "MINT_PYRAMID",        //!< LINEAR_PYRAMID
};

} /* namespace cell */
} /* namespace mint */
} /* namespace axom */

#endif /* CELLTYPE_HXX_ */
