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

typedef conduit::DataType::TypeID TypeID;



typedef int IDType;

const IDType InvalidID = -1;


template< int TYPEID >
TypeID getTypeID()
{
  return static_cast<TypeID>(TYPEID);
}



} /* end namespace sidre */
} /* end namespace asctoolkit */


#endif /* SIDRETYPES_HPP_ */
