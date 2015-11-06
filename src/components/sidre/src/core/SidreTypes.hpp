/**
 *  \file SidreTypes.hpp
 *
 *  \brief File containing types used in SiDRe toolkit component.
 *
 */

#ifndef SIDRETYPES_HPP_
#define SIDRETYPES_HPP_

// Other CS Toolkit headers
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
 * \brief TypeID is used to identify the type of a buffer (ATK_INT8_T, etc).
 */
typedef int TypeID;

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
inline bool isNameValid(const std::string& name)
{
  return name != InvalidName;
}



/*
 *************************************************************************
 *
 * Given a Sidre type enum create a Conduit DataType.
 *
 *************************************************************************
 */
#if 0
inline DataType createDataType( TypeID type, long len )
{
  DataType rval;

  switch( type )
  {
  case DataType::INT8_T:
    rval = DataType::int8(len);
    break;
  case DataType::INT16_T:
    rval = DataType::int16(len);
    break;
  case DataType::INT32_T:
    rval = DataType::int32(len);
    break;
  case DataType::INT64_T:
    rval = DataType::int64(len);
    break;
  case DataType::UINT8_T:
    rval = DataType::uint8(len);
    break;
  case DataType::UINT16_T:
    rval = DataType::uint16(len);
    break;
  case DataType::UINT32_T:
    rval = DataType::uint32(len);
    break;
  case DataType::UINT64_T:
    rval = DataType::uint64(len);
    break;
  case DataType::FLOAT32_T:
    rval = DataType::float32(len);
    break;
  case DataType::FLOAT64_T:
    rval = DataType::float64(len);
    break;
#if 0
  case DataType::CHAR8_STR_T:
    rval = DataType::c_char(len);
    break;
#endif
  default:
    break;
//      ATK_ERROR( "getTypeID(int) passed invalid type" );


  }

  return rval;

}
#endif



} /* end namespace sidre */
} /* end namespace asctoolkit */


#endif /* SIDRETYPES_HPP_ */
