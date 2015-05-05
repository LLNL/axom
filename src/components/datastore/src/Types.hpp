/**
 *  \file Types.hpp
 *
 *  \brief a file that contains some typedefs and type based functionality
 *
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <iostream>
#include "conduit/conduit.h"


namespace sidre
{
 
// The IDType should always be an unsigned type.  Our internal checks assume the var is >= 0.   
typedef conduit::index_t IDType;

#ifdef USE_CXX11

#else

#ifndef nullptr
#define nullptr 0
#endif

#endif

// typedef short int int16;
// typedef unsigned short int uint16;
// typedef long int int32;
// typedef unsigned long int uint32;
// typedef long long int int64;
// typedef unsigned long long int uint64;
// typedef std::size_t sizet;
//
// typedef float real32;
// typedef double real64;
//
//
//
//
// /**
//  * \class rtTypes
//  * \brief class that holds some basic type functionality
//  */
// namespace rtTypes
// {
//
//   /**
//    * \enum TypeID
//    * \brief enums to tag each basic type
//    */
//   enum TypeID
//   {
//     int16_id, //!< int16_id
//     uint16_id,//!< uint16_id
//     int32_id, //!< int32_id
//     uint32_id,//!< uint32_id
//     int64_id, //!< int64_id
//     uint64_id,//!< uint64_id
//     sizet_id, //!< sizet_id
//     real32_id,//!< real32_id
//     real64_id,//!< real64_id
//     undefined
//   };
//
//
//   template< typename TYPE >
//   inline TypeID GetTypeID()
//   {
//     return undefined;
//   }
//   template<> inline TypeID GetTypeID<int>() { return int32_id; }
//   template<> inline TypeID GetTypeID<int16>() { return int16_id; }
//   template<> inline TypeID GetTypeID<uint16>() { return uint16_id; }
//   template<> inline TypeID GetTypeID<int32>() { return int32_id; }
//   template<> inline TypeID GetTypeID<uint32>() { return uint32_id; }
//   template<> inline TypeID GetTypeID<int64>() { return int64_id; }
//   template<> inline TypeID GetTypeID<uint64>() { return uint64_id; }
// //  template<> inline TypeID GetTypeID<sizet>()  { return sizet_id; }
//
//   template<>  inline TypeID GetTypeID<real32>() { return real32_id; }
//   template<> inline TypeID GetTypeID<real64>() { return real64_id; }
//
//
//   struct sizeof_wrapper
//   {
//     template< typename T >
//     std::size_t get()
//     {
//       return sizeof(T);
//     }
//   } ;
//
//   inline std::size_t sizeofType( const TypeID typeID )
//   {
//     std::size_t rval;
//     switch( typeID )
//     {
//       case( int16_id ):
//       {
//         rval = sizeof(int16);
//         break;
//       }
//       case( uint16_id ):
//       {
//         rval = sizeof(uint16);
//         break;
//       }
//       case( int32_id ):
//       {
//         rval = sizeof(int32);
//         break;
//       }
//       case( uint32_id ):
//       {
//         rval = sizeof(uint32);
//         break;
//       }
//       case( int64_id ):
//       {
//         rval = sizeof(int64);
//         break;
//       }
//       case( uint64_id ):
//       {
//         rval = sizeof(uint64);
//         break;
//       }
//       case( sizet_id ):
//       {
//         rval = sizeof(sizet);
//         break;
//       }
//       case( real32_id ):
//       {
//         rval = sizeof(real32);
//         break;
//       }
//       case( real64_id ):
//       {
//         rval = sizeof(real64);
//         break;
//       }
//       default:
//       {
//         throw std::exception();
//       }
//     }
//
//     return rval;
//   }
// }
//
//
//
// /**
//  * struct that should encapsulate all information needed to describe a chunk of data.
//  */
// struct DataShape
// {
//   DataShape():
//     m_dataPtr(NULL),
//     m_numDimensions( 0 ),
//     m_dimensions(NULL),
//     m_dataStride(0)
//   {}
//
//   DataShape( const int& numDimensions, const std::size_t* const dimensions ):
//     m_dataPtr(NULL),
//     m_numDimensions( numDimensions ),
//     m_dimensions(NULL),
//     m_dataStride(0)
//   {
//     m_dimensions = new std::size_t[m_numDimensions];
//
//     for( int i=0 ; i<m_numDimensions ; ++i )
//     {
//       m_dimensions[i] = dimensions[i];
//     }
//     m_dataStride = m_dimensions[m_numDimensions-1];
//
//   }
//
//   DataShape( const int& length ):
//     m_dataPtr(NULL),
//     m_numDimensions( 1 ),
//     m_dimensions(NULL),
//     m_dataStride(0)
//   {
//     m_dimensions = new std::size_t[m_numDimensions];
//     m_dimensions[0] = length;
//     m_dataStride = m_dimensions[m_numDimensions-1];
//   }
//
//
//   DataShape( const DataShape& source):
//     m_dataPtr(source.m_dataPtr),
//     m_numDimensions( source.m_numDimensions ),
//     m_dimensions(NULL),
//     m_dataStride(source.m_dataStride)
//   {
//     m_dimensions = new std::size_t[m_numDimensions];
//     for( int i=0 ; i<m_numDimensions ; ++i )
//     {
//       m_dimensions[i] = source.m_dimensions[i];
//     }
//   }
//
//   /*
//   DataShape( DataShape&& source):
//     m_dataPtr(std::move(source.m_dataPtr)),
//     m_numDimensions( std::move(source.m_numDimensions) ),
//     m_dimensions(std::move(source.m_dimensions)),
//     m_dataStride(std::move(source.m_dataStride))
//   {
//     source.m_dataPtr = NULL;
//     source.m_numDimensions = 0;
//     source.m_dimensions = NULL;
//     source.m_dataStride = 0;
//
//   }
//   */
//
//   DataShape& operator=( const DataShape& source )
//   {
//     m_dataPtr = NULL;
//     m_numDimensions = source.m_numDimensions ;
//     m_dataStride = source.m_dataStride;
//
//     if( m_dimensions!=NULL )
//     {
//       delete[] m_dimensions;
//     }
//     m_dimensions = new std::size_t[m_numDimensions];
//
//     for( int i=0 ; i<m_numDimensions ; ++i )
//     {
//       m_dimensions[i] = source.m_dimensions[i];
//     }
//
//     return *this;
//   }
//
//   ~DataShape()
//   {
//     if( m_dimensions!=NULL )
//     {
//       delete[] m_dimensions;
//     }
//   }
//
//   void* m_dataPtr;
//   int m_numDimensions;
//   std::size_t* m_dimensions;
//   std::size_t m_dataStride;
//
// };



}


#endif /* TYPES_HPP_ */
