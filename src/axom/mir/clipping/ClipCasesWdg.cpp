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

int numClipCasesWdg = 64;

int numClipShapesWdg[64] = {
  1,  3,  3,  9,  3,  9,  9,  2, // cases 0 - 7
  3,  2,  15,  13,  15,  13,  10,  9, // cases 8 - 15
  3,  15,  2,  13,  15,  10,  13,  9, // cases 16 - 23
  9,  13,  13,  2,  10,  5,  5,  3, // cases 24 - 31
  3,  15,  15,  10,  2,  13,  13,  9, // cases 32 - 39
  9,  13,  10,  5,  13,  2,  5,  3, // cases 40 - 47
  9,  10,  13,  5,  13,  5,  2,  3, // cases 48 - 55
  2,  9,  9,  3,  9,  3,  3,  1  // cases 56 - 63
};

int startClipShapesWdg[64] = {
  0, 8, 29, 50, 115, 136, 201, 266, // cases 0 - 7
  282, 303, 321, 421, 508, 608, 695, 768, // cases 8 - 15
  833, 854, 954, 972, 1059, 1159, 1232, 1319, // cases 16 - 23
  1384, 1449, 1536, 1623, 1641, 1714, 1748, 1782, // cases 24 - 31
  1803, 1824, 1924, 2024, 2097, 2115, 2202, 2289, // cases 32 - 39
  2354, 2419, 2506, 2579, 2613, 2700, 2718, 2752, // cases 40 - 47
  2773, 2838, 2911, 2998, 3032, 3119, 3153, 3171, // cases 48 - 55
  3192, 3208, 3273, 3338, 3359, 3424, 3445, 3466  // cases 56 - 63
};

unsigned char clipShapesWdg[] = {
 // Case #0: Unique case #1
  ST_WDG, COLOR0, P0, P1, P2, P3, P4, P5, 
 // Case #1: Unique case #2
  ST_WDG, COLOR0, EA, EC, EG, P1, P2, P3, 
  ST_PYR, COLOR0, P1, P2, P5, P4, P3, 
  ST_TET, COLOR1, EG, EA, EC, P0, 
 // Case #2: (cloned #1)
  ST_WDG, COLOR0, EB, EA, EH, P2, P0, P4, 
  ST_PYR, COLOR0, P2, P0, P3, P5, P4, 
  ST_TET, COLOR1, EH, EB, EA, P1, 
 // Case #3: Unique case #3
  ST_PNT, 0, COLOR0, 7, P2, P3, P4, EB, EC, EG, EH, 
  ST_TET, COLOR0, P4, P5, P3, P2, 
  ST_TET, COLOR0, P2, P3, P4, N0, 
  ST_PYR, COLOR0, EG, EH, P4, P3, N0, 
  ST_PYR, COLOR0, EB, EH, EG, EC, N0, 
  ST_TET, COLOR0, P2, EB, EC, N0, 
  ST_PYR, COLOR0, P2, EC, EG, P3, N0, 
  ST_PYR, COLOR0, EH, EB, P2, P4, N0, 
  ST_WDG, COLOR1, EC, EG, P0, EB, EH, P1, 
 // Case #4: (cloned #1)
  ST_WDG, COLOR0, EC, EB, EI, P0, P1, P5, 
  ST_PYR, COLOR0, P0, P1, P4, P3, P5, 
  ST_TET, COLOR1, EI, EC, EB, P2, 
 // Case #5: (cloned #3)
  ST_PNT, 0, COLOR0, 7, P1, P5, P3, EA, EB, EI, EG, 
  ST_TET, COLOR0, P3, P4, P5, P1, 
  ST_TET, COLOR0, P1, P5, P3, N0, 
  ST_PYR, COLOR0, EI, EG, P3, P5, N0, 
  ST_PYR, COLOR0, EA, EG, EI, EB, N0, 
  ST_TET, COLOR0, P1, EA, EB, N0, 
  ST_PYR, COLOR0, P1, EB, EI, P5, N0, 
  ST_PYR, COLOR0, EG, EA, P1, P3, N0, 
  ST_WDG, COLOR1, EB, EI, P2, EA, EG, P0, 
 // Case #6: (cloned #3)
  ST_PNT, 0, COLOR0, 7, P0, P4, P5, EC, EA, EH, EI, 
  ST_TET, COLOR0, P5, P3, P4, P0, 
  ST_TET, COLOR0, P0, P4, P5, N0, 
  ST_PYR, COLOR0, EH, EI, P5, P4, N0, 
  ST_PYR, COLOR0, EC, EI, EH, EA, N0, 
  ST_TET, COLOR0, P0, EC, EA, N0, 
  ST_PYR, COLOR0, P0, EA, EH, P4, N0, 
  ST_PYR, COLOR0, EI, EC, P0, P5, N0, 
  ST_WDG, COLOR1, EA, EH, P1, EC, EI, P2, 
 // Case #7: Unique case #4
  ST_WDG, COLOR0, EG, EH, EI, P3, P4, P5, 
  ST_WDG, COLOR1, P0, P1, P2, EG, EH, EI, 
 // Case #8: (cloned #1)
  ST_WDG, COLOR0, P4, P5, P0, ED, EF, EG, 
  ST_PYR, COLOR0, P4, P1, P2, P5, P0, 
  ST_TET, COLOR1, EG, EF, ED, P3, 
 // Case #9: Unique case #5
  ST_HEX, COLOR0, P1, P2, P5, P4, EA, EC, EF, ED, 
  ST_WDG, COLOR1, P0, EA, EC, P3, ED, EF, 
 // Case #10: Unique case #6
  ST_PNT, 0, NOCOLOR, 6, EA, EB, EH, ED, EF, EG, 
  ST_PYR, COLOR0, P5, P0, EG, EF, N0, 
  ST_TET, COLOR0, P0, EA, EG, N0, 
  ST_PYR, COLOR0, P0, P2, EB, EA, N0, 
  ST_TET, COLOR0, P5, P2, P0, N0, 
  ST_PYR, COLOR0, P4, EH, EB, P2, N0, 
  ST_TET, COLOR0, P5, P4, P2, N0, 
  ST_PYR, COLOR0, EF, ED, P4, P5, N0, 
  ST_TET, COLOR0, ED, EH, P4, N0, 
  ST_PYR, COLOR1, EG, EA, P1, P3, N0, 
  ST_PYR, COLOR1, P3, P1, EH, ED, N0, 
  ST_TET, COLOR1, P3, ED, EF, N0, 
  ST_TET, COLOR1, EF, EG, P3, N0, 
  ST_TET, COLOR1, P1, EB, EH, N0, 
  ST_TET, COLOR1, P1, EA, EB, N0, 
 // Case #11: Unique case #7
  ST_PNT, 0, NOCOLOR, 5, EB, EC, EF, ED, EH, 
  ST_PYR, COLOR0, P4, P5, EF, ED, N0, 
  ST_TET, COLOR0, ED, EH, P4, N0, 
  ST_PYR, COLOR0, EC, EF, P5, P2, N0, 
  ST_PYR, COLOR0, EB, P2, P4, EH, N0, 
  ST_TET, COLOR0, P4, P2, P5, N0, 
  ST_TET, COLOR0, P2, EB, EC, N0, 
  ST_TET, COLOR1, P0, P1, P3, N0, 
  ST_PYR, COLOR1, EC, EB, P1, P0, N0, 
  ST_PYR, COLOR1, EC, P0, P3, EF, N0, 
  ST_TET, COLOR1, EF, P3, ED, N0, 
  ST_PYR, COLOR1, P3, P1, EH, ED, N0, 
  ST_TET, COLOR1, P1, EB, EH, N0, 
 // Case #12: (cloned #10)
  ST_PNT, 0, NOCOLOR, 6, EF, ED, EG, EC, EB, EI, 
  ST_PYR, COLOR0, P1, EB, EI, P5, N0, 
  ST_TET, COLOR0, P5, EI, EF, N0, 
  ST_PYR, COLOR0, P5, EF, ED, P4, N0, 
  ST_TET, COLOR0, P1, P5, P4, N0, 
  ST_PYR, COLOR0, P0, P4, ED, EG, N0, 
  ST_TET, COLOR0, P1, P4, P0, N0, 
  ST_PYR, COLOR0, EB, P1, P0, EC, N0, 
  ST_TET, COLOR0, EC, P0, EG, N0, 
  ST_PYR, COLOR1, EI, P2, P3, EF, N0, 
  ST_PYR, COLOR1, P2, EC, EG, P3, N0, 
  ST_TET, COLOR1, P2, EB, EC, N0, 
  ST_TET, COLOR1, EB, P2, EI, N0, 
  ST_TET, COLOR1, P3, EG, ED, N0, 
  ST_TET, COLOR1, P3, ED, EF, N0, 
 // Case #13: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EB, EA, ED, EF, EI, 
  ST_PYR, COLOR0, P5, EF, ED, P4, N0, 
  ST_TET, COLOR0, EF, P5, EI, N0, 
  ST_PYR, COLOR0, EA, P1, P4, ED, N0, 
  ST_PYR, COLOR0, EB, EI, P5, P1, N0, 
  ST_TET, COLOR0, P5, P4, P1, N0, 
  ST_TET, COLOR0, P1, EA, EB, N0, 
  ST_TET, COLOR1, P0, P3, P2, N0, 
  ST_PYR, COLOR1, EA, P0, P2, EB, N0, 
  ST_PYR, COLOR1, EA, ED, P3, P0, N0, 
  ST_TET, COLOR1, ED, EF, P3, N0, 
  ST_PYR, COLOR1, P3, EF, EI, P2, N0, 
  ST_TET, COLOR1, P2, EI, EB, N0, 
 // Case #14: Unique case #8
  ST_PNT, 0, COLOR1, 7, ED, EF, EI, EH, P3, P2, P1, 
  ST_TET, COLOR0, P0, EC, EA, EG, 
  ST_WDG, COLOR0, EF, EI, P5, ED, EH, P4, 
  ST_WDG, COLOR1, P2, P1, P3, EC, EA, EG, 
  ST_PYR, COLOR1, EF, ED, EH, EI, N0, 
  ST_PYR, COLOR1, EH, P1, P2, EI, N0, 
  ST_TET, COLOR1, P3, P2, P1, N0, 
  ST_TET, COLOR1, P3, ED, EF, N0, 
  ST_PYR, COLOR1, ED, P3, P1, EH, N0, 
  ST_PYR, COLOR1, EI, P2, P3, EF, N0, 
 // Case #15: Unique case #9
  ST_PNT, 0, COLOR1, 7, P1, P2, P3, EF, ED, EH, EI, 
  ST_WDG, COLOR0, ED, P4, EH, EF, P5, EI, 
  ST_TET, COLOR1, P0, P2, P1, P3, 
  ST_PYR, COLOR1, EF, ED, EH, EI, N0, 
  ST_PYR, COLOR1, EI, EH, P1, P2, N0, 
  ST_TET, COLOR1, P3, P2, P1, N0, 
  ST_TET, COLOR1, P3, ED, EF, N0, 
  ST_PYR, COLOR1, P3, P1, EH, ED, N0, 
  ST_PYR, COLOR1, P2, P3, EF, EI, N0, 
 // Case #16: (cloned #1)
  ST_WDG, COLOR0, P5, P3, P1, EE, ED, EH, 
  ST_PYR, COLOR0, P5, P2, P0, P3, P1, 
  ST_TET, COLOR1, EH, ED, EE, P4, 
 // Case #17: (cloned #10)
  ST_PNT, 0, NOCOLOR, 6, ED, EE, EH, EA, EC, EG, 
  ST_PYR, COLOR0, P2, EC, EG, P3, N0, 
  ST_TET, COLOR0, P3, EG, ED, N0, 
  ST_PYR, COLOR0, P3, ED, EE, P5, N0, 
  ST_TET, COLOR0, P2, P3, P5, N0, 
  ST_PYR, COLOR0, P1, P5, EE, EH, N0, 
  ST_TET, COLOR0, P2, P5, P1, N0, 
  ST_PYR, COLOR0, EC, P2, P1, EA, N0, 
  ST_TET, COLOR0, EA, P1, EH, N0, 
  ST_PYR, COLOR1, EG, P0, P4, ED, N0, 
  ST_PYR, COLOR1, P0, EA, EH, P4, N0, 
  ST_TET, COLOR1, P0, EC, EA, N0, 
  ST_TET, COLOR1, EC, P0, EG, N0, 
  ST_TET, COLOR1, P4, EH, EE, N0, 
  ST_TET, COLOR1, P4, EE, ED, N0, 
 // Case #18: (cloned #9)
  ST_HEX, COLOR0, P2, P0, P3, P5, EB, EA, ED, EE, 
  ST_WDG, COLOR1, P1, EB, EA, P4, EE, ED, 
 // Case #19: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EC, EB, EE, ED, EG, 
  ST_PYR, COLOR0, P3, ED, EE, P5, N0, 
  ST_TET, COLOR0, ED, P3, EG, N0, 
  ST_PYR, COLOR0, EB, P2, P5, EE, N0, 
  ST_PYR, COLOR0, EC, EG, P3, P2, N0, 
  ST_TET, COLOR0, P3, P5, P2, N0, 
  ST_TET, COLOR0, P2, EB, EC, N0, 
  ST_TET, COLOR1, P1, P4, P0, N0, 
  ST_PYR, COLOR1, EB, P1, P0, EC, N0, 
  ST_PYR, COLOR1, EB, EE, P4, P1, N0, 
  ST_TET, COLOR1, EE, ED, P4, N0, 
  ST_PYR, COLOR1, P4, ED, EG, P0, N0, 
  ST_TET, COLOR1, P0, EG, EC, N0, 
 // Case #20: (cloned #10)
  ST_PNT, 0, NOCOLOR, 6, EB, EC, EI, EE, ED, EH, 
  ST_PYR, COLOR0, P3, P1, EH, ED, N0, 
  ST_TET, COLOR0, P1, EB, EH, N0, 
  ST_PYR, COLOR0, P1, P0, EC, EB, N0, 
  ST_TET, COLOR0, P3, P0, P1, N0, 
  ST_PYR, COLOR0, P5, EI, EC, P0, N0, 
  ST_TET, COLOR0, P3, P5, P0, N0, 
  ST_PYR, COLOR0, ED, EE, P5, P3, N0, 
  ST_TET, COLOR0, EE, EI, P5, N0, 
  ST_PYR, COLOR1, EH, EB, P2, P4, N0, 
  ST_PYR, COLOR1, P4, P2, EI, EE, N0, 
  ST_TET, COLOR1, P4, EE, ED, N0, 
  ST_TET, COLOR1, ED, EH, P4, N0, 
  ST_TET, COLOR1, P2, EC, EI, N0, 
  ST_TET, COLOR1, P2, EB, EC, N0, 
 // Case #21: (cloned #14)
  ST_PNT, 0, COLOR1, 7, EE, ED, EG, EI, P4, P0, P2, 
  ST_TET, COLOR0, P1, EA, EB, EH, 
  ST_WDG, COLOR0, ED, EG, P3, EE, EI, P5, 
  ST_WDG, COLOR1, P0, P2, P4, EA, EB, EH, 
  ST_PYR, COLOR1, ED, EE, EI, EG, N0, 
  ST_PYR, COLOR1, EI, P2, P0, EG, N0, 
  ST_TET, COLOR1, P4, P0, P2, N0, 
  ST_TET, COLOR1, P4, EE, ED, N0, 
  ST_PYR, COLOR1, EE, P4, P2, EI, N0, 
  ST_PYR, COLOR1, EG, P0, P4, ED, N0, 
 // Case #22: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EC, EA, ED, EE, EI, 
  ST_PYR, COLOR0, P5, P3, ED, EE, N0, 
  ST_TET, COLOR0, EE, EI, P5, N0, 
  ST_PYR, COLOR0, EA, ED, P3, P0, N0, 
  ST_PYR, COLOR0, EC, P0, P5, EI, N0, 
  ST_TET, COLOR0, P5, P0, P3, N0, 
  ST_TET, COLOR0, P0, EC, EA, N0, 
  ST_TET, COLOR1, P1, P2, P4, N0, 
  ST_PYR, COLOR1, EA, EC, P2, P1, N0, 
  ST_PYR, COLOR1, EA, P1, P4, ED, N0, 
  ST_TET, COLOR1, ED, P4, EE, N0, 
  ST_PYR, COLOR1, P4, P2, EI, EE, N0, 
  ST_TET, COLOR1, P2, EC, EI, N0, 
 // Case #23: (cloned #15)
  ST_PNT, 0, COLOR1, 7, P2, P0, P4, ED, EE, EI, EG, 
  ST_WDG, COLOR0, EE, P5, EI, ED, P3, EG, 
  ST_TET, COLOR1, P1, P0, P2, P4, 
  ST_PYR, COLOR1, ED, EE, EI, EG, N0, 
  ST_PYR, COLOR1, EG, EI, P2, P0, N0, 
  ST_TET, COLOR1, P4, P0, P2, N0, 
  ST_TET, COLOR1, P4, EE, ED, N0, 
  ST_PYR, COLOR1, P4, P2, EI, EE, N0, 
  ST_PYR, COLOR1, P0, P4, ED, EG, N0, 
 // Case #24: (cloned #3)
  ST_PNT, 0, COLOR0, 7, P5, P0, P1, EE, EF, EG, EH, 
  ST_TET, COLOR0, P1, P0, P2, P5, 
  ST_TET, COLOR0, P5, P1, P0, N0, 
  ST_PYR, COLOR0, EG, P0, P1, EH, N0, 
  ST_PYR, COLOR0, EE, EF, EG, EH, N0, 
  ST_TET, COLOR0, P5, EF, EE, N0, 
  ST_PYR, COLOR0, P5, P0, EG, EF, N0, 
  ST_PYR, COLOR0, EH, P1, P5, EE, N0, 
  ST_WDG, COLOR1, EE, EH, P4, EF, EG, P3, 
 // Case #25: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EE, EF, EC, EA, EH, 
  ST_PYR, COLOR0, P1, EA, EC, P2, N0, 
  ST_TET, COLOR0, EA, P1, EH, N0, 
  ST_PYR, COLOR0, EF, P5, P2, EC, N0, 
  ST_PYR, COLOR0, EE, EH, P1, P5, N0, 
  ST_TET, COLOR0, P1, P2, P5, N0, 
  ST_TET, COLOR0, P5, EF, EE, N0, 
  ST_TET, COLOR1, P3, P0, P4, N0, 
  ST_PYR, COLOR1, EF, P3, P4, EE, N0, 
  ST_PYR, COLOR1, EF, EC, P0, P3, N0, 
  ST_TET, COLOR1, EC, EA, P0, N0, 
  ST_PYR, COLOR1, P0, EA, EH, P4, N0, 
  ST_TET, COLOR1, P4, EH, EE, N0, 
 // Case #26: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EF, EE, EB, EA, EG, 
  ST_PYR, COLOR0, P0, P2, EB, EA, N0, 
  ST_TET, COLOR0, EA, EG, P0, N0, 
  ST_PYR, COLOR0, EE, EB, P2, P5, N0, 
  ST_PYR, COLOR0, EF, P5, P0, EG, N0, 
  ST_TET, COLOR0, P0, P5, P2, N0, 
  ST_TET, COLOR0, P5, EF, EE, N0, 
  ST_TET, COLOR1, P4, P3, P1, N0, 
  ST_PYR, COLOR1, EE, EF, P3, P4, N0, 
  ST_PYR, COLOR1, EE, P4, P1, EB, N0, 
  ST_TET, COLOR1, EB, P1, EA, N0, 
  ST_PYR, COLOR1, P1, P3, EG, EA, N0, 
  ST_TET, COLOR1, P3, EF, EG, N0, 
 // Case #27: Unique case #10
  ST_WDG, COLOR0, EF, P5, EE, EC, P2, EB, 
  ST_HEX, COLOR1, P3, P4, EE, EF, P0, P1, EB, EC, 
 // Case #28: (cloned #14)
  ST_PNT, 0, COLOR1, 7, EC, EB, EH, EG, P2, P4, P3, 
  ST_TET, COLOR0, P5, EF, EE, EI, 
  ST_WDG, COLOR0, EC, EG, P0, EB, EH, P1, 
  ST_WDG, COLOR1, EE, EF, EI, P4, P3, P2, 
  ST_PYR, COLOR1, EB, EH, EG, EC, N0, 
  ST_PYR, COLOR1, EG, EH, P4, P3, N0, 
  ST_TET, COLOR1, P2, P3, P4, N0, 
  ST_TET, COLOR1, P2, EB, EC, N0, 
  ST_PYR, COLOR1, EC, EG, P3, P2, N0, 
  ST_PYR, COLOR1, EH, EB, P2, P4, N0, 
 // Case #29: Unique case #11
  ST_TET, COLOR0, P1, EA, EB, EH, 
  ST_TET, COLOR0, EF, EE, P5, EI, 
  ST_WDG, COLOR1, P2, P3, P4, EI, EF, EE, 
  ST_TET, COLOR1, P2, P3, P4, P0, 
  ST_WDG, COLOR1, P2, P4, P0, EB, EH, EA, 
 // Case #30: (cloned #29)
  ST_TET, COLOR0, P5, EF, EE, EI, 
  ST_TET, COLOR0, EA, P0, EC, EG, 
  ST_WDG, COLOR1, EG, EA, EC, P3, P1, P2, 
  ST_TET, COLOR1, P3, P2, P1, P4, 
  ST_WDG, COLOR1, EF, EI, EE, P3, P2, P4, 
 // Case #31: Unique case #12
  ST_TET, COLOR0, EF, EI, EE, P5, 
  ST_WDG, COLOR1, EI, EE, EF, P2, P4, P3, 
  ST_PYR, COLOR1, P0, P1, P4, P3, P2, 
 // Case #32: (cloned #1)
  ST_WDG, COLOR0, P3, P4, P2, EF, EE, EI, 
  ST_PYR, COLOR0, P3, P0, P1, P4, P2, 
  ST_TET, COLOR1, EI, EE, EF, P5, 
 // Case #33: (cloned #10)
  ST_PNT, 0, NOCOLOR, 6, EC, EA, EG, EF, EE, EI, 
  ST_PYR, COLOR0, P4, P2, EI, EE, N0, 
  ST_TET, COLOR0, P2, EC, EI, N0, 
  ST_PYR, COLOR0, P2, P1, EA, EC, N0, 
  ST_TET, COLOR0, P4, P1, P2, N0, 
  ST_PYR, COLOR0, P3, EG, EA, P1, N0, 
  ST_TET, COLOR0, P4, P3, P1, N0, 
  ST_PYR, COLOR0, EE, EF, P3, P4, N0, 
  ST_TET, COLOR0, EF, EG, P3, N0, 
  ST_PYR, COLOR1, EI, EC, P0, P5, N0, 
  ST_PYR, COLOR1, P5, P0, EG, EF, N0, 
  ST_TET, COLOR1, P5, EF, EE, N0, 
  ST_TET, COLOR1, EE, EI, P5, N0, 
  ST_TET, COLOR1, P0, EA, EG, N0, 
  ST_TET, COLOR1, P0, EC, EA, N0, 
 // Case #34: (cloned #10)
  ST_PNT, 0, NOCOLOR, 6, EE, EF, EI, EB, EA, EH, 
  ST_PYR, COLOR0, P0, EA, EH, P4, N0, 
  ST_TET, COLOR0, P4, EH, EE, N0, 
  ST_PYR, COLOR0, P4, EE, EF, P3, N0, 
  ST_TET, COLOR0, P0, P4, P3, N0, 
  ST_PYR, COLOR0, P2, P3, EF, EI, N0, 
  ST_TET, COLOR0, P0, P3, P2, N0, 
  ST_PYR, COLOR0, EA, P0, P2, EB, N0, 
  ST_TET, COLOR0, EB, P2, EI, N0, 
  ST_PYR, COLOR1, EH, P1, P5, EE, N0, 
  ST_PYR, COLOR1, P1, EB, EI, P5, N0, 
  ST_TET, COLOR1, P1, EA, EB, N0, 
  ST_TET, COLOR1, EA, P1, EH, N0, 
  ST_TET, COLOR1, P5, EI, EF, N0, 
  ST_TET, COLOR1, P5, EF, EE, N0, 
 // Case #35: (cloned #14)
  ST_PNT, 0, COLOR1, 7, EF, EE, EH, EG, P5, P1, P0, 
  ST_TET, COLOR0, P2, EB, EC, EI, 
  ST_WDG, COLOR0, EE, EH, P4, EF, EG, P3, 
  ST_WDG, COLOR1, P1, P0, P5, EB, EC, EI, 
  ST_PYR, COLOR1, EE, EF, EG, EH, N0, 
  ST_PYR, COLOR1, EG, P0, P1, EH, N0, 
  ST_TET, COLOR1, P5, P1, P0, N0, 
  ST_TET, COLOR1, P5, EF, EE, N0, 
  ST_PYR, COLOR1, EF, P5, P0, EG, N0, 
  ST_PYR, COLOR1, EH, P1, P5, EE, N0, 
 // Case #36: (cloned #9)
  ST_HEX, COLOR0, P0, P1, P4, P3, EC, EB, EE, EF, 
  ST_WDG, COLOR1, P2, EC, EB, P5, EF, EE, 
 // Case #37: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EA, EB, EE, EF, EG, 
  ST_PYR, COLOR0, P3, P4, EE, EF, N0, 
  ST_TET, COLOR0, EF, EG, P3, N0, 
  ST_PYR, COLOR0, EB, EE, P4, P1, N0, 
  ST_PYR, COLOR0, EA, P1, P3, EG, N0, 
  ST_TET, COLOR0, P3, P1, P4, N0, 
  ST_TET, COLOR0, P1, EA, EB, N0, 
  ST_TET, COLOR1, P2, P0, P5, N0, 
  ST_PYR, COLOR1, EB, EA, P0, P2, N0, 
  ST_PYR, COLOR1, EB, P2, P5, EE, N0, 
  ST_TET, COLOR1, EE, P5, EF, N0, 
  ST_PYR, COLOR1, P5, P0, EG, EF, N0, 
  ST_TET, COLOR1, P0, EA, EG, N0, 
 // Case #38: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EA, EC, EF, EE, EH, 
  ST_PYR, COLOR0, P4, EE, EF, P3, N0, 
  ST_TET, COLOR0, EE, P4, EH, N0, 
  ST_PYR, COLOR0, EC, P0, P3, EF, N0, 
  ST_PYR, COLOR0, EA, EH, P4, P0, N0, 
  ST_TET, COLOR0, P4, P3, P0, N0, 
  ST_TET, COLOR0, P0, EC, EA, N0, 
  ST_TET, COLOR1, P2, P5, P1, N0, 
  ST_PYR, COLOR1, EC, P2, P1, EA, N0, 
  ST_PYR, COLOR1, EC, EF, P5, P2, N0, 
  ST_TET, COLOR1, EF, EE, P5, N0, 
  ST_PYR, COLOR1, P5, EE, EH, P1, N0, 
  ST_TET, COLOR1, P1, EH, EA, N0, 
 // Case #39: (cloned #15)
  ST_PNT, 0, COLOR1, 7, P0, P1, P5, EE, EF, EG, EH, 
  ST_WDG, COLOR0, EF, P3, EG, EE, P4, EH, 
  ST_TET, COLOR1, P2, P1, P0, P5, 
  ST_PYR, COLOR1, EE, EF, EG, EH, N0, 
  ST_PYR, COLOR1, EH, EG, P0, P1, N0, 
  ST_TET, COLOR1, P5, P1, P0, N0, 
  ST_TET, COLOR1, P5, EF, EE, N0, 
  ST_PYR, COLOR1, P5, P0, EG, EF, N0, 
  ST_PYR, COLOR1, P1, P5, EE, EH, N0, 
 // Case #40: (cloned #3)
  ST_PNT, 0, COLOR0, 7, P4, P2, P0, ED, EE, EI, EG, 
  ST_TET, COLOR0, P0, P2, P1, P4, 
  ST_TET, COLOR0, P4, P0, P2, N0, 
  ST_PYR, COLOR0, EI, P2, P0, EG, N0, 
  ST_PYR, COLOR0, ED, EE, EI, EG, N0, 
  ST_TET, COLOR0, P4, EE, ED, N0, 
  ST_PYR, COLOR0, P4, P2, EI, EE, N0, 
  ST_PYR, COLOR0, EG, P0, P4, ED, N0, 
  ST_WDG, COLOR1, ED, EG, P3, EE, EI, P5, 
 // Case #41: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EE, ED, EA, EC, EI, 
  ST_PYR, COLOR0, P2, P1, EA, EC, N0, 
  ST_TET, COLOR0, EC, EI, P2, N0, 
  ST_PYR, COLOR0, ED, EA, P1, P4, N0, 
  ST_PYR, COLOR0, EE, P4, P2, EI, N0, 
  ST_TET, COLOR0, P2, P4, P1, N0, 
  ST_TET, COLOR0, P4, EE, ED, N0, 
  ST_TET, COLOR1, P3, P5, P0, N0, 
  ST_PYR, COLOR1, ED, EE, P5, P3, N0, 
  ST_PYR, COLOR1, ED, P3, P0, EA, N0, 
  ST_TET, COLOR1, EA, P0, EC, N0, 
  ST_PYR, COLOR1, P0, P5, EI, EC, N0, 
  ST_TET, COLOR1, P5, EE, EI, N0, 
 // Case #42: (cloned #14)
  ST_PNT, 0, COLOR1, 7, EB, EA, EG, EI, P1, P3, P5, 
  ST_TET, COLOR0, P4, EE, ED, EH, 
  ST_WDG, COLOR0, EB, EI, P2, EA, EG, P0, 
  ST_WDG, COLOR1, ED, EE, EH, P3, P5, P1, 
  ST_PYR, COLOR1, EA, EG, EI, EB, N0, 
  ST_PYR, COLOR1, EI, EG, P3, P5, N0, 
  ST_TET, COLOR1, P1, P5, P3, N0, 
  ST_TET, COLOR1, P1, EA, EB, N0, 
  ST_PYR, COLOR1, EB, EI, P5, P1, N0, 
  ST_PYR, COLOR1, EG, EA, P1, P3, N0, 
 // Case #43: (cloned #29)
  ST_TET, COLOR0, P4, EE, ED, EH, 
  ST_TET, COLOR0, EC, P2, EB, EI, 
  ST_WDG, COLOR1, EI, EC, EB, P5, P0, P1, 
  ST_TET, COLOR1, P5, P1, P0, P3, 
  ST_WDG, COLOR1, EE, EH, ED, P5, P1, P3, 
 // Case #44: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, ED, EE, EB, EC, EG, 
  ST_PYR, COLOR0, P0, EC, EB, P1, N0, 
  ST_TET, COLOR0, EC, P0, EG, N0, 
  ST_PYR, COLOR0, EE, P4, P1, EB, N0, 
  ST_PYR, COLOR0, ED, EG, P0, P4, N0, 
  ST_TET, COLOR0, P0, P1, P4, N0, 
  ST_TET, COLOR0, P4, EE, ED, N0, 
  ST_TET, COLOR1, P5, P2, P3, N0, 
  ST_PYR, COLOR1, EE, P5, P3, ED, N0, 
  ST_PYR, COLOR1, EE, EB, P2, P5, N0, 
  ST_TET, COLOR1, EB, EC, P2, N0, 
  ST_PYR, COLOR1, P2, EC, EG, P3, N0, 
  ST_TET, COLOR1, P3, EG, ED, N0, 
 // Case #45: (cloned #27)
  ST_WDG, COLOR0, EE, P4, ED, EB, P1, EA, 
  ST_HEX, COLOR1, P5, P3, ED, EE, P2, P0, EA, EB, 
 // Case #46: (cloned #29)
  ST_TET, COLOR0, P0, EC, EA, EG, 
  ST_TET, COLOR0, EE, ED, P4, EH, 
  ST_WDG, COLOR1, P1, P5, P3, EH, EE, ED, 
  ST_TET, COLOR1, P1, P5, P3, P2, 
  ST_WDG, COLOR1, P1, P3, P2, EA, EG, EC, 
 // Case #47: (cloned #31)
  ST_TET, COLOR0, EE, EH, ED, P4, 
  ST_WDG, COLOR1, EH, ED, EE, P1, P3, P5, 
  ST_PYR, COLOR1, P2, P0, P3, P5, P1, 
 // Case #48: (cloned #3)
  ST_PNT, 0, COLOR0, 7, P3, P1, P2, EF, ED, EH, EI, 
  ST_TET, COLOR0, P2, P1, P0, P3, 
  ST_TET, COLOR0, P3, P2, P1, N0, 
  ST_PYR, COLOR0, EH, P1, P2, EI, N0, 
  ST_PYR, COLOR0, EF, ED, EH, EI, N0, 
  ST_TET, COLOR0, P3, ED, EF, N0, 
  ST_PYR, COLOR0, P3, P1, EH, ED, N0, 
  ST_PYR, COLOR0, EI, P2, P3, EF, N0, 
  ST_WDG, COLOR1, EF, EI, P5, ED, EH, P4, 
 // Case #49: (cloned #14)
  ST_PNT, 0, COLOR1, 7, EA, EC, EI, EH, P0, P5, P4, 
  ST_TET, COLOR0, P3, ED, EF, EG, 
  ST_WDG, COLOR0, EA, EH, P1, EC, EI, P2, 
  ST_WDG, COLOR1, EF, ED, EG, P5, P4, P0, 
  ST_PYR, COLOR1, EC, EI, EH, EA, N0, 
  ST_PYR, COLOR1, EH, EI, P5, P4, N0, 
  ST_TET, COLOR1, P0, P4, P5, N0, 
  ST_TET, COLOR1, P0, EC, EA, N0, 
  ST_PYR, COLOR1, EA, EH, P4, P0, N0, 
  ST_PYR, COLOR1, EI, EC, P0, P5, N0, 
 // Case #50: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, EF, ED, EA, EB, EI, 
  ST_PYR, COLOR0, P2, EB, EA, P0, N0, 
  ST_TET, COLOR0, EB, P2, EI, N0, 
  ST_PYR, COLOR0, ED, P3, P0, EA, N0, 
  ST_PYR, COLOR0, EF, EI, P2, P3, N0, 
  ST_TET, COLOR0, P2, P0, P3, N0, 
  ST_TET, COLOR0, P3, ED, EF, N0, 
  ST_TET, COLOR1, P4, P1, P5, N0, 
  ST_PYR, COLOR1, ED, P4, P5, EF, N0, 
  ST_PYR, COLOR1, ED, EA, P1, P4, N0, 
  ST_TET, COLOR1, EA, EB, P1, N0, 
  ST_PYR, COLOR1, P1, EB, EI, P5, N0, 
  ST_TET, COLOR1, P5, EI, EF, N0, 
 // Case #51: (cloned #29)
  ST_TET, COLOR0, P2, EB, EC, EI, 
  ST_TET, COLOR0, ED, EF, P3, EG, 
  ST_WDG, COLOR1, P0, P4, P5, EG, ED, EF, 
  ST_TET, COLOR1, P0, P4, P5, P1, 
  ST_WDG, COLOR1, P0, P5, P1, EC, EI, EB, 
 // Case #52: (cloned #11)
  ST_PNT, 0, NOCOLOR, 5, ED, EF, EC, EB, EH, 
  ST_PYR, COLOR0, P1, P0, EC, EB, N0, 
  ST_TET, COLOR0, EB, EH, P1, N0, 
  ST_PYR, COLOR0, EF, EC, P0, P3, N0, 
  ST_PYR, COLOR0, ED, P3, P1, EH, N0, 
  ST_TET, COLOR0, P1, P3, P0, N0, 
  ST_TET, COLOR0, P3, ED, EF, N0, 
  ST_TET, COLOR1, P5, P4, P2, N0, 
  ST_PYR, COLOR1, EF, ED, P4, P5, N0, 
  ST_PYR, COLOR1, EF, P5, P2, EC, N0, 
  ST_TET, COLOR1, EC, P2, EB, N0, 
  ST_PYR, COLOR1, P2, P4, EH, EB, N0, 
  ST_TET, COLOR1, P4, ED, EH, N0, 
 // Case #53: (cloned #29)
  ST_TET, COLOR0, P3, ED, EF, EG, 
  ST_TET, COLOR0, EB, P1, EA, EH, 
  ST_WDG, COLOR1, EH, EB, EA, P4, P2, P0, 
  ST_TET, COLOR1, P4, P0, P2, P5, 
  ST_WDG, COLOR1, ED, EG, EF, P4, P0, P5, 
 // Case #54: (cloned #27)
  ST_WDG, COLOR0, ED, P3, EF, EA, P0, EC, 
  ST_HEX, COLOR1, P4, P5, EF, ED, P1, P2, EC, EA, 
 // Case #55: (cloned #31)
  ST_TET, COLOR0, ED, EG, EF, P3, 
  ST_WDG, COLOR1, EG, EF, ED, P0, P5, P4, 
  ST_PYR, COLOR1, P1, P2, P5, P4, P0, 
 // Case #56: (cloned #7)
  ST_WDG, COLOR0, P0, P1, P2, EG, EH, EI, 
  ST_WDG, COLOR1, EG, EH, EI, P3, P4, P5, 
 // Case #57: (cloned #15)
  ST_PNT, 0, COLOR1, 7, P4, P5, P0, EC, EA, EH, EI, 
  ST_WDG, COLOR0, EC, P2, EI, EA, P1, EH, 
  ST_TET, COLOR1, P3, P4, P5, P0, 
  ST_PYR, COLOR1, EC, EI, EH, EA, N0, 
  ST_PYR, COLOR1, EI, P5, P4, EH, N0, 
  ST_TET, COLOR1, P0, P4, P5, N0, 
  ST_TET, COLOR1, P0, EC, EA, N0, 
  ST_PYR, COLOR1, P0, EA, EH, P4, N0, 
  ST_PYR, COLOR1, P5, EI, EC, P0, N0, 
 // Case #58: (cloned #15)
  ST_PNT, 0, COLOR1, 7, P5, P3, P1, EA, EB, EI, EG, 
  ST_WDG, COLOR0, EA, P0, EG, EB, P2, EI, 
  ST_TET, COLOR1, P4, P5, P3, P1, 
  ST_PYR, COLOR1, EA, EG, EI, EB, N0, 
  ST_PYR, COLOR1, EG, P3, P5, EI, N0, 
  ST_TET, COLOR1, P1, P5, P3, N0, 
  ST_TET, COLOR1, P1, EA, EB, N0, 
  ST_PYR, COLOR1, P1, EB, EI, P5, N0, 
  ST_PYR, COLOR1, P3, EG, EA, P1, N0, 
 // Case #59: (cloned #31)
  ST_TET, COLOR0, EC, EB, EI, P2, 
  ST_WDG, COLOR1, P5, P1, P0, EI, EB, EC, 
  ST_PYR, COLOR1, P3, P0, P1, P4, P5, 
 // Case #60: (cloned #15)
  ST_PNT, 0, COLOR1, 7, P3, P4, P2, EB, EC, EG, EH, 
  ST_WDG, COLOR0, EB, P1, EH, EC, P0, EG, 
  ST_TET, COLOR1, P5, P3, P4, P2, 
  ST_PYR, COLOR1, EB, EH, EG, EC, N0, 
  ST_PYR, COLOR1, EH, P4, P3, EG, N0, 
  ST_TET, COLOR1, P2, P3, P4, N0, 
  ST_TET, COLOR1, P2, EB, EC, N0, 
  ST_PYR, COLOR1, P2, EC, EG, P3, N0, 
  ST_PYR, COLOR1, P4, EH, EB, P2, N0, 
 // Case #61: (cloned #31)
  ST_TET, COLOR0, EB, EA, EH, P1, 
  ST_WDG, COLOR1, P4, P0, P2, EH, EA, EB, 
  ST_PYR, COLOR1, P5, P2, P0, P3, P4, 
 // Case #62: (cloned #31)
  ST_TET, COLOR0, EA, EC, EG, P0, 
  ST_WDG, COLOR1, P3, P2, P1, EG, EC, EA, 
  ST_PYR, COLOR1, P4, P1, P2, P5, P3, 
 // Case #63: Unique case #13
  ST_WDG, COLOR1, P0, P1, P2, P3, P4, P5, 
 // Dummy
  0
};

//---------------------------------------------------------------------------
// Axom modifications
const size_t clipShapesWdgSize = sizeof(clipShapesWdg) / sizeof(unsigned char);
} // namespace visit
} // namespace clipping
} // namespace mir
} // namespace axom
// clang-format on
//---------------------------------------------------------------------------
