/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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
 *  \file SidreTypes.hpp
 *
 *  \brief File containing types used in the Sidre component of axom.
 *
 */

#ifndef SIDRETYPES_HPP_
#define SIDRETYPES_HPP_

#include "SidreDataTypeIds.h"
#include "conduit.hpp"
#include "axom/core/Types.hpp"

namespace axom
{
namespace sidre
{
// Some typedefs to make Conduit usage easier and less visible in the Sidre API.

/*!
 * \brief The Node class is the primary object in Conduit.
 */
typedef conduit::Node Node;

/*!
 * \brief DataType is a general Conduit descriptor.
 */
typedef conduit::DataType DataType;

/*!
 * \brief A Conduit Schema describes the data in a Node
 */
typedef conduit::Schema Schema;


/*!
 * \brief IndexType is used for any labeling of a sidre object by an
 *        integer identifier.
 */
typedef SIDRE_IndexType IndexType;

/*!
 * \brief Common invalid index identifier used in sidre.
 */
const IndexType InvalidIndex = SIDRE_InvalidIndex;

/*!
 * \brief Returns true if idx is valid, else false.
 *
 * Used for the loop test when iterating over the
 * Buffers or Attributes in a DataStore or the Views or Groups in a Group.
 */
inline bool indexIsValid(IndexType idx)
{
  return idx != InvalidIndex;
}

/*!
 * \brief Common invalid name (string) identifier used in sidre.
 */
const std::string InvalidName;

/*!
 * \brief Returns true if name is valid, else false.
 */
inline bool nameIsValid(const std::string& name)
{
  return name != InvalidName;
}

/*!
 * \brief Enum that holds the numeric data type id options for sidre types.
 */

enum DataTypeId
{
  NO_TYPE_ID = SIDRE_NO_TYPE_ID,
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

/// @cond INCLUDE_DETAIL

/*!
 * \brief The detail namespace contains code that is either used internally by
 *  the sidre implementation or is under evaluation.
 */
namespace detail
{

/*!
 * \brief Type traits to assist in converting compiler types to the appropriate
 *  data type ids.
 */
template<typename T> struct SidreTT
{
  static const DataTypeId id = NO_TYPE_ID;
};

template<> struct SidreTT<common::int8>
{
  static const DataTypeId id = INT8_ID;
};
template<> struct SidreTT<common::int16>
{
  static const DataTypeId id = INT16_ID;
};
template<> struct SidreTT<common::int32>
{
  static const DataTypeId id = INT32_ID;
};
template<> struct SidreTT<common::int64>
{
  static const DataTypeId id = INT64_ID;
};

template<> struct SidreTT<common::uint8>
{
  static const DataTypeId id = UINT8_ID;
};
template<> struct SidreTT<common::uint16>
{
  static const DataTypeId id = UINT16_ID;
};
template<> struct SidreTT<common::uint32>
{
  static const DataTypeId id = UINT32_ID;
};
template<> struct SidreTT<common::uint64>
{
  static const DataTypeId id = UINT64_ID;
};

template<> struct SidreTT<common::float32>
{
  static const DataTypeId id = FLOAT32_ID;
};
template<> struct SidreTT<common::float64>
{
  static const DataTypeId id = FLOAT64_ID;
};
} /* end namespace detail */
/// @endcond

/*!
 * \brief TypeID is used to identify the type of a buffer (SIDRE_INT8_ID, etc).
 */
typedef DataTypeId TypeID;

/*!
 * \brief Convenience function to convert int to TypeID type.
 *
 *  Used to convert C defines to C++ enumerations.
 */
inline TypeID getTypeID( const int typeID )
{
  return static_cast<TypeID>(typeID);
}

} /* end namespace sidre */
} /* end namespace axom */


#endif /* SIDRETYPES_HPP_ */
