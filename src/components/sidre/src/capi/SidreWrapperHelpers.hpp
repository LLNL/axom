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
 *  \file SidreWrapperHelpers.hpp
 *
 *  \brief File used to contain helper functions for C API wrappers.
 *         User code should not include this file.
 *
 */

#ifndef SIDREWRAPPERHELPERS_HPP_
#define SIDREWRAPPERHELPERS_HPP_

#include "conduit/Bitwidth_Style_Types.h"

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
  TypeID rval = CONDUIT_EMPTY_T;

  switch( typeID )
  {
#if 0
  case 0:
    rval = CONDUIT_EMPTY_T;
    break;
  case 1:
    rval = CONDUIT_OBJECT_T;
    break;
  case 2:
    rval = CONDUIT_LIST_T;
    break;
#endif
  case ATK_INT8_T:
    rval = CONDUIT_INT8_T;
    break;
  case ATK_INT16_T:
    rval = CONDUIT_INT16_T;
    break;
  case ATK_INT32_T:
    rval = CONDUIT_INT32_T;
    break;
  case ATK_INT64_T:
    rval = CONDUIT_INT64_T;
    break;
  case ATK_UINT8_T:
    rval = CONDUIT_UINT8_T;
    break;
  case ATK_UINT16_T:
    rval = CONDUIT_UINT16_T;
    break;
  case ATK_UINT32_T:
    rval = CONDUIT_UINT32_T;
    break;
  case ATK_UINT64_T:
    rval = CONDUIT_UINT64_T;
    break;
  case ATK_FLOAT32_T:
    rval = CONDUIT_FLOAT32_T;
    break;
  case ATK_FLOAT64_T:
    rval = CONDUIT_FLOAT64_T;
    break;
  case ATK_CHAR8_STR_T:
    rval = CONDUIT_CHAR8_STR_T;
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
    rval = CONDUIT_EMPTY_T;
    break;
//      ATK_ERROR( "getTypeID(int) passed invalid type" );


  }

  return rval;

}


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* SIDREWRAPPERHELPERS_HPP_ */



