/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


//
// DataTypes.h
//
// Data enums used by C and Fortran
#ifndef SIDRE_DATATYPES_H
#define SIDRE_DATATYPES_H

typedef int ATK_IndexType;

typedef long int ATK_SidreLength;

#define ATK_InvalidIndex  -1

#define ATK_InvalidName   NULL
//
//  Must be kept in sync with sidretypes.inc
#define ATK_INT8_T   3
#define ATK_INT16_T   4
#define ATK_INT32_T   5
#define ATK_INT64_T   6
#define ATK_UINT8_T   7
#define ATK_UINT16_T   8
#define ATK_UINT32_T   9
#define ATK_UINT64_T   10
#define ATK_FLOAT32_T   11
#define ATK_FLOAT64_T    12
#define ATK_CHAR8_STR_T   13

#define ATK_C_INT_T   14
#define ATK_C_LONG_T   15
#define ATK_C_FLOAT_T   16
#define ATK_C_DOUBLE_T   17


#endif  // SIDRE_DATATYPES_H
