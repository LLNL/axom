/**
 *  \file CommonTypes.hpp
 *
 *  \brief File containing common types used in toolkit components.
 *
 */

#ifndef COMMONTYPES_HPP_
#define COMMONTYPES_HPP_

namespace asctoolkit
{
 
namespace common
{

#ifdef USE_CXX11
#define ATK_NULLPTR nullptr
#else
#define ATK_NULLPTR NULL
#endif


} /* end namespace common */
} /* end namespace asctoolkit */


#endif /* COMMONTYPES_HPP_ */
