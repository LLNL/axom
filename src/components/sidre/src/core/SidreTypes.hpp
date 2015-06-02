/**
 *  \file SidreTypes.hpp
 *
 *  \brief File containing types used in SiDRe toolkit component.
 *
 */

#ifndef SIDRETYPES_HPP_
#define SIDRETYPES_HPP_

// Other CS Toolkit headers
#include "conduit/conduit.hpp"

// SiDRe project headers
#include "sidre/SidreTypes.h"

namespace asctoolkit
{
namespace sidre
{

typedef conduit::DataType DataType;
typedef conduit::DataType::TypeID TypeID;


/*!
 * \brief IndexType is used for any labeling of a sidre object by an
 *        integer identifier.
 */
typedef int IndexType;

/*!
 * \brief Common invalid index identifier used in sidre.
 */
const IndexType InvalidIndex = -1;


/*!
 * \brief Common invalid name (string) identifier used in sidre.
 */
const std::string InvalidName;
///
inline bool isNameValid(const std::string& name)
{
  return name != InvalidName;
}

template< int TYPEID >
inline TypeID getTypeID()
{
  return static_cast<TypeID>(TYPEID);
}

inline TypeID getTypeID( const ATK_TypeEnum typeID )
{
  TypeID rval = DataType::EMPTY_T;

  switch( typeID )
  {
#if 0
  case 0:
    rval = DataType::EMPTY_T;
    break;
  case 1:
    rval = DataType::OBJECT_T;
    break;
  case 2:
    rval = DataType::LIST_T;
    break;
#endif
  case ATK_INT8_T:
    rval = DataType::INT8_T;
    break;
  case ATK_INT16_T:
    rval = DataType::INT16_T;
    break;
  case ATK_INT32_T:
    rval = DataType::INT32_T;
    break;
  case ATK_INT64_T:
    rval = DataType::INT64_T;
    break;
  case ATK_UINT8_T:
    rval = DataType::UINT8_T;
    break;
  case ATK_UINT16_T:
    rval = DataType::UINT16_T;
    break;
  case ATK_UINT32_T:
    rval = DataType::UINT32_T;
    break;
  case ATK_UINT64_T:
    rval = DataType::UINT64_T;
    break;
  case ATK_FLOAT32_T:
    rval = DataType::FLOAT32_T;
    break;
  case ATK_FLOAT64_T:
    rval = DataType::FLOAT64_T;
    break;
  case ATK_CHAR8_STR_T:
    rval = DataType::CHAR8_STR_T;
    break;
  default:
    rval = DataType::EMPTY_T;
    break;
//      ATK_ERROR( "getTypeID(int) passed invalid type" );


  }

  return rval;

}

/*
 *************************************************************************
 *
 * Given a Sidre type enum create a Conduit DataType.
 *
 *************************************************************************
 */
inline conduit::DataType createConduitDataType( const TypeID type, long len )
{
  conduit::DataType rval;

  switch( type )
  {
  case DataType::INT8_T:
    rval = conduit::DataType::int8(len);
    break;
  case DataType::INT16_T:
    rval = conduit::DataType::int16(len);
    break;
  case DataType::INT32_T:
    rval = conduit::DataType::int32(len);
    break;
  case DataType::INT64_T:
    rval = conduit::DataType::int64(len);
    break;
  case DataType::UINT8_T:
    rval = conduit::DataType::uint8(len);
    break;
  case DataType::UINT16_T:
    rval = conduit::DataType::uint16(len);
    break;
  case DataType::UINT32_T:
    rval = conduit::DataType::uint32(len);
    break;
  case DataType::UINT64_T:
    rval = conduit::DataType::uint64(len);
    break;
  case DataType::FLOAT32_T:
    rval = conduit::DataType::float32(len);
    break;
  case DataType::FLOAT64_T:
    rval = conduit::DataType::float64(len);
    break;
#if 0
  case DataType::CHAR8_STR_T:
    rval = conduit::DataType::c_char(len);
    break;
#endif
  default:
    break;
//      ATK_ERROR( "getTypeID(int) passed invalid type" );


  }

  return rval;

}



} /* end namespace sidre */
} /* end namespace asctoolkit */


#endif /* SIDRETYPES_HPP_ */
