/**
 *  \file Types.hpp
 *
 *  \brief a file that contains some typedefs and type based functionality
 *
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

// TODO:  We should pull out the types out of conduit and put them here...
// -- Aaron
#include "boost/unordered_map.hpp"
#include "conduit/conduit.h"

namespace asctoolkit
{
 
namespace common
{

// The IDType should always be an unsigned type.  Our internal checks assume the var is >= 0.
// //TODO: We should define our own IDType instead of using conduit's. -- Aaron
typedef conduit::index_t IDType;

// Add typedefs for C++11 only container types.  We are supporting using C++11 containers if
// there is an equivalent boost library solution we can use on compilers that don't support C++11.
#ifdef USE_CXX11
typedef boost::unordered_map< std::string, common::IDType> UnorderedMapStringToIDType;
#else
typedef boost::unordered_map< std::string, common::IDType> UnorderedMapStringToIDType;
#endif

#ifdef USE_CXX11
#define ATK_NULLPTR nullptr
#else

const // It is a const object...
class atk_nullptr_t 
{
  public:
    template<class T>
    inline operator T*() const // convertible to any type of null non-member pointer...
    { return 0; }
 
    template<class C, class T>
    inline operator T C::*() const   // or any type of null member pointer...
    { return 0; }
 
  private:
    void operator&() const;  // Can't take address of nullptr
 
} atk_nullptr = {};

#define ATK_NULLPTR atk_nullptr

#endif


} /* end namespace common */
} /* end namespace asctoolkit */


#endif /* TYPES_HPP_ */
