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



} /* end namespace sidre */
} /* end namespace asctoolkit */


#endif /* SIDRETYPES_HPP_ */
