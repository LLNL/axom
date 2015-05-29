/**
 *  \file SidreTypes.hpp
 *
 *  \brief File containing types used in SiDRe toolkit component.
 *
 */

#ifndef SIDRETYPES_HPP_
#define SIDRETYPES_HPP_

namespace asctoolkit
{
namespace sidre
{

typedef conduit::DataType DataType;
typedef conduit::DataType::TypeID TypeID;



typedef int IDType;

const IDType InvalidID = -1;

#include "sidre/DataTypes.h"


template< int TYPEID >
inline TypeID getTypeID()
{
  return static_cast<TypeID>(TYPEID);
}

inline TypeID getTypeID( const int typeID )
{
  TypeID rval = DataType::EMPTY_T;

  switch( typeID )
  {
    case 0:
      rval = DataType::EMPTY_T;
      break;
    case 1:
      rval = DataType::OBJECT_T;
      break;
    case 2:
      rval = DataType::LIST_T;
      break;
    case 3:
      rval = DataType::INT8_T;
      break;
    case 4:
      rval = DataType::INT16_T;
      break;
    case 5:
      rval = DataType::INT32_T;
      break;
    case 6:
      rval = DataType::INT64_T;
      break;
    case 7:
      rval = DataType::UINT8_T;
      break;
    case 8:
      rval = DataType::UINT16_T;
      break;
    case 9:
      rval = DataType::UINT32_T;
      break;
    case 10:
      rval = DataType::UINT64_T;
      break;
    case 11:
      rval = DataType::FLOAT32_T;
      break;
    case 12:
      rval = DataType::FLOAT64_T;
      break;
    case 13:
      rval = DataType::CHAR8_STR_T;
      break;
    default:
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
inline conduit::DataType createConduitDataType( const ATK_TypeID type, long len )
{
 conduit::DataType rval;

  switch( type )
  {
  case INT8_T:
      rval = conduit::DataType::int8(len);
      break;
  case INT16_T:
      rval = conduit::DataType::int16(len);
      break;
  case INT32_T:
      rval = conduit::DataType::int32(len);
      break;
  case INT64_T:
      rval = conduit::DataType::int64(len);
      break;
  case UINT8_T:
      rval = conduit::DataType::uint8(len);
      break;
  case UINT16_T:
      rval = conduit::DataType::uint16(len);
      break;
  case UINT32_T:
      rval = conduit::DataType::uint32(len);
      break;
  case UINT64_T:
      rval = conduit::DataType::uint64(len);
      break;
  case FLOAT32_T:
      rval = conduit::DataType::float32(len);
      break;
  case FLOAT64_T:
      rval = conduit::DataType::float64(len);
      break;
#if 0
  case CHAR8_STR_T:
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
