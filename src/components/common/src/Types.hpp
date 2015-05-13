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
// //TODO: Should we should define our own IDType instead of using conduit's? -- Aaron
typedef conduit::index_t IDType;

#ifdef USE_CXX11
#define ATK_NULLPTR nullptr
#else
// Note:
// When casting this type to another pointer type, you will need to explicitly cast it with a static_cast.
// The compiler will throw an error if any implicit void* casts are detected.  This allows us to prevent
// implicit casts due to developer error.
#define ATK_NULLPTR (void*)0
#endif


} /* end namespace common */
} /* end namespace asctoolkit */


#endif /* TYPES_HPP_ */
