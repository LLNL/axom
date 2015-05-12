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
#include "conduit/conduit.h"

namespace asctoolkit
{
 
namespace common
{

// The IDType should always be an unsigned type.  Our internal checks assume the var is >= 0.
// //TODO: We should define our own IDType instead of using conduit's. -- Aaron
typedef conduit::index_t IDType;

#ifdef USE_CXX11
#define ATK_NULLPTR nullptr
#else
//#define ATK_NULLPTR (void*)0

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

#define ATK_NULLPTR asctoolkit::common::atk_nullptr

#endif


} /* end namespace common */
} /* end namespace asctoolkit */


#endif /* TYPES_HPP_ */
