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
#define ATK_NULLPTR (void*)0
#endif


} /* end namespace common */
} /* end namespace asctoolkit */


#endif /* TYPES_HPP_ */
