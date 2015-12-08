/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file CellType.hxx
 *
 * \date Sep 12, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef CELLTYPE_HXX_
#define CELLTYPE_HXX_

namespace meshtk {


enum CellType {
  UNDEFINED = -1,

  VERTEX,
  LINE,

  LINEAR_TRIANGLE,
  LINEAR_QUAD,
  LINEAR_TET,
  LINEAR_HEX,
  LINEAR_PRISM,
  LINEAR_PYRAMID,

  MIXED,
  NUM_CELL_TYPES
};

namespace cell {

static const int vtk_types[] = {
  1,    // VERTEX          -> VTK_VERTEX
  3,    // LINE            -> VTK_LINE
  5,    // LINEAR_TRIANGLE -> VTK_TRIANGLE
  9,    // LINEAR_QUAD     -> VTK_QUAD
  10,   // LINEAR_TET      -> VTK_TET
  12,   // LINEAR_HEX      -> VTK_HEXAHEDRON
  13,   // LINEAR_PRISM    -> VTK_WEDGE
  14,   // LINEAR_PYRAMID  -> VTK_PYRAMID
  0     // MIXED           -> VTK_EMPTY_CELL
};

static const int num_nodes[] = {
   1,   // VERTEX
   2,   // LINE
   3,   // LINEAR_TRIANGLE
   4,   // LINEAR_QUAD
   4,   // LINEAR_TET
   8,   // LINEAR_HEX
   6,   // LINEAR_PRISM
   5,   // LINEAR_PYRAMID

   5    // MIXED
};

} /* namespace cell */

} /* namespace meshtk */

#endif /* CELLTYPE_HXX_ */
