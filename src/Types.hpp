/**
 *  \file Types.hpp
 *
 *  \brief a file that contains some typedefs and type based functionality
 *
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_
#include <typeindex>
#include <iostream>
/**
 *
 */
namespace DataStore
{
typedef short int int16;
typedef unsigned short int uint16;
typedef long int int32;
typedef unsigned long int uint32;
typedef long long int int64;
typedef unsigned long long int uint64;
typedef std::size_t sizet;

typedef float real32;
typedef double real64;





/**
 * \class rtTypes
 * \brief class that holds some basic type functionality
 */
namespace rtTypes
{

  /**
   * \enum TypeID
   * \brief enums to tag each basic type
   */
  enum class TypeID
  {
    int16_id, //!< int16_id
    uint16_id,//!< uint16_id
    int32_id, //!< int32_id
    uint32_id,//!< uint32_id
    int64_id, //!< int64_id
    uint64_id,//!< uint64_id
    sizet_id, //!< sizet_id
    real32_id,//!< real32_id
    real64_id,//!< real64_id
    undefined
  };


  template< typename TYPE >
  inline TypeID GetTypeID()
  {
    return TypeID::undefined;
  }
  template<> inline TypeID GetTypeID<int>() { return TypeID::int32_id; }
  template<> inline TypeID GetTypeID<int16>() { return TypeID::int16_id; }
  template<> inline TypeID GetTypeID<uint16>() { return TypeID::uint16_id; }
  template<> inline TypeID GetTypeID<int32>() { return TypeID::int32_id; }
  template<> inline TypeID GetTypeID<uint32>() { return TypeID::uint32_id; }
  template<> inline TypeID GetTypeID<int64>() { return TypeID::int64_id; }
  template<> inline TypeID GetTypeID<uint64>() { return TypeID::uint64_id; }
//  template<> inline TypeID GetTypeID<sizet>()  { return TypeID::sizet_id; }

  template<>  inline TypeID GetTypeID<real32>() { return TypeID::real32_id; }
  template<> inline TypeID GetTypeID<real64>() { return TypeID::real64_id; }



  template< typename TYPE >
  inline TYPE* CastPtr( void* const ptr )
  {
    return static_cast<TYPE*>(ptr);
  }





  /**
   * @tparam WRAPPER a wrapper type for the function that you want to call. WRAPPER must have member get() that calls the
   *                 actual function that is to be called.
   * @tparam ArgsF variadic argument typelist for the function to be called
   * @param type[in]  what TypeID
   * @param wrapper[in] the wrapper
   * @param args[in/out]  the variadic arguments
   * @return auto depends on the return type of wrapper.get()
   *
   * \brief This function will call
   */
  template< typename WRAPPER, typename...ArgsF >
  inline auto SelectFunction( const TypeID type, WRAPPER& wrapper, ArgsF&... args ) -> decltype( wrapper.template get<double>(args...) )
  {
    decltype(wrapper.template get<double>(args...)) rval;
    switch( type )
    {
      case( TypeID::int16_id ):
      {
        rval = wrapper.template get<int16>(args... );
        break;
      }
      case( TypeID::uint16_id ):
      {
        rval = wrapper.template get<uint16>(args... );
        break;
      }
      case( TypeID::int32_id ):
      {
        rval = wrapper.template get<int32>(args... );
        break;
      }
      case( TypeID::uint32_id ):
      {
        rval = wrapper.template get<uint32>(args... );
        break;
      }
      case( TypeID::int64_id ):
      {
        rval = wrapper.template get<int64>(args... );
        break;
      }
      case( TypeID::uint64_id ):
      {
        rval = wrapper.template get<uint64>(args... );
        break;
      }
      case( TypeID::sizet_id ):
      {
        rval = wrapper.template get<sizet>(args... );
        break;
      }
      case( TypeID::real32_id ):
      {
        rval =  wrapper.template get<real32>(args...);
        break;
      }
      case( TypeID::real64_id ):
      {
        rval =  wrapper.template get<real64>(args...);
        break;
      }
      default:
      {
        throw std::exception();
      }
    }

    return rval;
  }


  struct sizeof_wrapper
  {
    template< typename T >
    std::size_t get()
    {
      return sizeof(T);
    }
  } ;

  inline std::size_t sizeofType( const TypeID typeID )
  {
    sizeof_wrapper sizeofWrapper;
    return SelectFunction( typeID, sizeofWrapper );
  }


}




/**
 * struct that should encapsulate all information needed to describe a chunk of data.
 */
struct DataShape
{
  DataShape():
    m_dataPtr(nullptr),
    m_numDimensions( 0 ),
    m_dimensions(nullptr),
    m_dataStride(0)
  {}

  DataShape( const int& numDimensions, const std::size_t* const dimensions ):
    m_dataPtr(nullptr),
    m_numDimensions( numDimensions ),
    m_dimensions(nullptr),
    m_dataStride(0)
  {
    m_dimensions = new std::size_t[m_numDimensions];

    for( int i=0 ; i<m_numDimensions ; ++i )
    {
      m_dimensions[i] = dimensions[i];
    }
    m_dataStride = m_dimensions[m_numDimensions-1];

  }

  DataShape( const int& length ):
    m_dataPtr(nullptr),
    m_numDimensions( 1 ),
    m_dimensions(nullptr),
    m_dataStride(0)
  {
    m_dimensions = new std::size_t[m_numDimensions];
    m_dimensions[0] = length;
    m_dataStride = m_dimensions[m_numDimensions-1];
  }


  DataShape( const DataShape& source):
    m_dataPtr(source.m_dataPtr),
    m_numDimensions( source.m_numDimensions ),
    m_dimensions(nullptr),
    m_dataStride(source.m_dataStride)
  {
    m_dimensions = new std::size_t[m_numDimensions];
    for( int i=0 ; i<m_numDimensions ; ++i )
    {
      m_dimensions[i] = source.m_dimensions[i];
    }
  }

  DataShape( DataShape&& source):
    m_dataPtr(std::move(source.m_dataPtr)),
    m_numDimensions( std::move(source.m_numDimensions) ),
    m_dimensions(std::move(source.m_dimensions)),
    m_dataStride(std::move(source.m_dataStride))
  {
    source.m_dataPtr = nullptr;
    source.m_numDimensions = 0;
    source.m_dimensions = nullptr;
    source.m_dataStride = 0;

  }


  ~DataShape()
  {
    if( m_dimensions!=nullptr )
    {
      delete[] m_dimensions;
    }
  }

  void* m_dataPtr;
  int m_numDimensions;
  std::size_t* m_dimensions;
  std::size_t m_dataStride;

};



}


#endif /* TYPES_HPP_ */
