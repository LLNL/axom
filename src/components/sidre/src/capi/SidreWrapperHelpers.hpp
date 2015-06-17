/**
 *  \file SidreWrapperHelpers.hpp
 *
 *  \brief File used to contain helper functions for C API wrappers.
 *         User code should not include this file.
 *
 */

#ifndef SIDREWRAPPERHELPERS_HPP_
#define SIDREWRAPPERHELPERS_HPP_


// SiDRe project headers
#include "sidre/SidreTypes.hpp"
#include "sidre/SidreTypes.h"

namespace asctoolkit
{
namespace sidre
{

template< int TYPEID >
inline TypeID getTypeID()
{
  return static_cast<TypeID>(TYPEID);
}

/*
 *************************************************************************
 *
 * Convert C #define value to TypeID.
 * Used with C wrappers.
 *
 *************************************************************************
 */
inline TypeID getTypeID( const int typeID )
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

  case ATK_C_INT_T:
    rval = CONDUIT_NATIVE_INT_DATATYPE_ID;
    break;
  case ATK_C_LONG_T:
    rval = CONDUIT_NATIVE_LONG_DATATYPE_ID;
    break;
  case ATK_C_FLOAT_T:
    rval = CONDUIT_NATIVE_FLOAT_DATATYPE_ID;
    break;
  case ATK_C_DOUBLE_T:
    rval = CONDUIT_NATIVE_DOUBLE_DATATYPE_ID;
    break;

  default:
    rval = DataType::EMPTY_T;
    break;
//      ATK_ERROR( "getTypeID(int) passed invalid type" );


  }

  return rval;

}


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* SIDREWRAPPERHELPERS_HPP_ */



