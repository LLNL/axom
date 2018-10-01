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


/**
 *  \file SiderDataTypeIds.h
 *
 *  \brief Contains defines used to establish type ids for numeric types in
 *   Sidre component.
 *
 */

#ifndef SIDREDATATYPEIDS_H_
#define SIDREDATATYPEIDS_H_

// Libraries and other axom headers
#include "conduit.h"

#include <stdint.h>   /* for int64_t */

typedef int64_t SIDRE_IndexType;

const SIDRE_IndexType SIDRE_InvalidIndex = -1;

#define SIDRE_NO_TYPE_ID    CONDUIT_EMPTY_ID
#define SIDRE_INT8_ID       CONDUIT_INT8_ID
#define SIDRE_INT16_ID      CONDUIT_INT16_ID
#define SIDRE_INT32_ID      CONDUIT_INT32_ID
#define SIDRE_INT64_ID      CONDUIT_INT64_ID
#define SIDRE_UINT8_ID      CONDUIT_UINT8_ID
#define SIDRE_UINT16_ID     CONDUIT_UINT16_ID
#define SIDRE_UINT32_ID     CONDUIT_UINT32_ID
#define SIDRE_UINT64_ID     CONDUIT_UINT64_ID
#define SIDRE_FLOAT32_ID    CONDUIT_FLOAT32_ID
#define SIDRE_FLOAT64_ID    CONDUIT_FLOAT64_ID
#define SIDRE_CHAR8_STR_ID  CONDUIT_CHAR8_STR_ID

#define SIDRE_INT_ID     CONDUIT_NATIVE_INT_ID
#define SIDRE_UINT_ID    CONDUIT_NATIVE_UNSIGNED_INT_ID
#define SIDRE_LONG_ID    CONDUIT_NATIVE_LONG_ID
#define SIDRE_ULONG_ID   CONDUIT_NATIVE_UNSIGNED_LONG_ID
#define SIDRE_FLOAT_ID   CONDUIT_NATIVE_FLOAT_ID
#define SIDRE_DOUBLE_ID  CONDUIT_NATIVE_DOUBLE_ID

#endif /* SIDREDATATYPEIDS_H_ */
