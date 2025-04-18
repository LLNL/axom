// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#include "ClipCases.h"
//---------------------------------------------------------------------------
// Axom modifications
// clang-format off
namespace axom {
namespace mir {
namespace clipping {
namespace visit {
//---------------------------------------------------------------------------

// Programmer: Jeremy Meredith
// Date      : August 11, 2003
//
// Modifications:
//    Jeremy Meredith, Mon Sep 15 17:30:21 PDT 2003
//    Added ability for Centroid-Points to have an associated color.

// This file is meant to be read and created by a program other than a
// compiler.  If you must modify it by hand, at least be nice to the 
// parser and don't add anything else to this file or rearrange it.

int numClipCasesTet = 16;

int numClipShapesTet[16] = {
  1,  2,  2,  2,  2,  2,  2,  2, // cases 0 - 7
  2,  2,  2,  2,  2,  2,  2,  1  // cases 8 - 15
};

int startClipShapesTet[16] = {
  0, 6, 20, 34, 50, 64, 80, 96, // cases 0 - 7
  110, 124, 140, 156, 170, 186, 200, 214  // cases 8 - 15
};

unsigned char clipShapesTet[] = {
 // Case #0: Unique case #1
  ST_TET, COLOR0, P0, P1, P2, P3, 
 // Case #1: Unique case #2
  ST_WDG, COLOR0, EA, ED, EC, P1, P3, P2, 
  ST_TET, COLOR1, P0, EA, EC, ED, 
 // Case #2: (cloned #1)
  ST_WDG, COLOR0, P0, P3, P2, EA, EE, EB, 
  ST_TET, COLOR1, P1, EB, EA, EE, 
 // Case #3: Unique case #3
  ST_WDG, COLOR0, P3, ED, EE, P2, EC, EB, 
  ST_WDG, COLOR1, P0, ED, EC, P1, EE, EB, 
 // Case #4: (cloned #1)
  ST_WDG, COLOR0, EC, EF, EB, P0, P3, P1, 
  ST_TET, COLOR1, P2, EC, EB, EF, 
 // Case #5: (cloned #3)
  ST_WDG, COLOR0, P1, EA, EB, P3, ED, EF, 
  ST_WDG, COLOR1, P2, EF, EB, P0, ED, EA, 
 // Case #6: (cloned #3)
  ST_WDG, COLOR0, P3, EE, EF, P0, EA, EC, 
  ST_WDG, COLOR1, P1, EE, EA, P2, EF, EC, 
 // Case #7: Unique case #4
  ST_TET, COLOR0, ED, EE, EF, P3, 
  ST_WDG, COLOR1, ED, EE, EF, P0, P1, P2, 
 // Case #8: (cloned #1)
  ST_WDG, COLOR0, P0, P2, P1, ED, EF, EE, 
  ST_TET, COLOR1, P3, EE, ED, EF, 
 // Case #9: (cloned #3)
  ST_WDG, COLOR0, P2, EC, EF, P1, EA, EE, 
  ST_WDG, COLOR1, P0, EC, EA, P3, EF, EE, 
 // Case #10: (cloned #3)
  ST_WDG, COLOR0, P0, EA, ED, P2, EB, EF, 
  ST_WDG, COLOR1, P3, EF, ED, P1, EB, EA, 
 // Case #11: (cloned #7)
  ST_TET, COLOR0, EC, EF, EB, P2, 
  ST_WDG, COLOR1, P0, P1, P3, EC, EB, EF, 
 // Case #12: (cloned #3)
  ST_WDG, COLOR0, P1, EB, EE, P0, EC, ED, 
  ST_WDG, COLOR1, P2, EB, EC, P3, EE, ED, 
 // Case #13: (cloned #7)
  ST_TET, COLOR0, EA, EB, EE, P1, 
  ST_WDG, COLOR1, EA, EB, EE, P0, P2, P3, 
 // Case #14: (cloned #7)
  ST_TET, COLOR0, EA, ED, EC, P0, 
  ST_WDG, COLOR1, P1, P2, P3, EA, EC, ED, 
 // Case #15: Unique case #5
  ST_TET, COLOR1, P0, P1, P2, P3, 
 // Dummy
  0
};

//---------------------------------------------------------------------------
// Axom modifications
const size_t clipShapesTetSize = sizeof(clipShapesTet) / sizeof(unsigned char);
} // namespace visit
} // namespace clipping
} // namespace mir
} // namespace axom
// clang-format on
//---------------------------------------------------------------------------
