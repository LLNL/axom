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
 *  \file SidreTypes.hpp
 *
 *  \brief File containing types used in SiDRe toolkit component.
 *
 */

#ifndef SIDRETYPES_HPP_
#define SIDRETYPES_HPP_

#include "SidreDataTypeIds.h"
#include "conduit.hpp"

namespace asctoolkit
{
namespace sidre
{

/*!
 * \brief DataType is a general Conduit descriptor.
 */
typedef conduit::DataType DataType;


/*!
 * \brief IndexType is used for any labeling of a sidre object by an
 *        integer identifier.
 */
typedef int IndexType;

typedef long int SidreLength;

/*!
 * \brief Common invalid index identifier used in sidre.
 */
const IndexType InvalidIndex = -1;
///
inline bool indexIsValid(IndexType idx)
{
  return idx != InvalidIndex;
}

/*!
 * \brief Common invalid name (string) identifier used in sidre.
 */
const std::string InvalidName;
///
inline bool nameIsValid(const std::string& name)
{
  return name != InvalidName;
}

/*!
 * \brief Enum that holds the numeric data type id options for sidre types.
 */

enum DataTypeId
{    
    EMPTY_ID = SIDRE_EMPTY_ID,
    INT8_ID  = SIDRE_INT8_ID,
    INT16_ID = SIDRE_INT16_ID,
    INT32_ID = SIDRE_INT32_ID,
    INT64_ID = SIDRE_INT64_ID,
    
    UINT8_ID  = SIDRE_UINT8_ID,
    UINT16_ID = SIDRE_UINT16_ID,
    UINT32_ID = SIDRE_UINT32_ID,
    UINT64_ID = SIDRE_UINT64_ID,
    
    FLOAT32_ID   = SIDRE_FLOAT32_ID,
    FLOAT64_ID   = SIDRE_FLOAT64_ID,
    CHAR8_STR_ID = SIDRE_CHAR8_STR_ID,
    
    INT_ID = SIDRE_INT_ID,
    UINT_ID = SIDRE_UINT_ID,
    LONG_ID = SIDRE_LONG_ID,
    ULONG_ID = SIDRE_ULONG_ID,
    FLOAT_ID = SIDRE_FLOAT_ID,
    DOUBLE_ID = SIDRE_DOUBLE_ID
};


/*!
 * \brief TypeID is used to identify the type of a buffer (SIDRE_INT8_ID, etc).
 */
typedef DataTypeId TypeID;



} /* end namespace sidre */
} /* end namespace asctoolkit */


#endif /* SIDRETYPES_HPP_ */
