/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


/**
 *  \file SiderDataTypeIds.h
 *
 *  \brief Contains defines used to establish type ids for numeric types in SiDRe toolkit component.
 *
 */

#ifndef SIDREDATATYPEIDS_H_
#define SIDREDATATYPEIDS_H_

// Other CS Toolkit headers
#include "conduit.hpp"

#define SIDRE_EMPTY_ID      CONDUIT_EMPTY_T
#define SIDRE_INT8_ID       CONDUIT_INT8_T
#define SIDRE_INT16_ID      CONDUIT_INT16_T
#define SIDRE_INT32_ID      CONDUIT_INT32_T
#define SIDRE_INT64_ID      CONDUIT_INT64_T
#define SIDRE_UINT8_ID      CONDUIT_UINT8_T
#define SIDRE_UINT16_ID     CONDUIT_UINT16_T
#define SIDRE_UINT32_ID     CONDUIT_UINT32_T
#define SIDRE_UINT64_ID     CONDUIT_UINT64_T
#define SIDRE_FLOAT32_ID    CONDUIT_FLOAT32_T
#define SIDRE_FLOAT64_ID    CONDUIT_FLOAT64_T
#define SIDRE_CHAR8_STR_ID  CONDUIT_CHAR8_STR_T

#define SIDRE_INT_ID     CONDUIT_NATIVE_INT_DATATYPE_ID
#define SIDRE_UINT_ID    CONDUIT_NATIVE_UNSIGNED_INT_DATATYPE_ID
#define SIDRE_LONG_ID    CONDUIT_NATIVE_LONG_DATATYPE_ID
#define SIDRE_ULONG_ID   CONDUIT_NATIVE_UNSIGNED_LONG_DATATYPE_ID
#define SIDRE_FLOAT_ID   CONDUIT_NATIVE_FLOAT_DATATYPE_ID
#define SIDRE_DOUBLE_ID  CONDUIT_NATIVE_DOUBLE_DATATYPE_ID

#endif /* SIDREDATATYPEIDS_H_ */

