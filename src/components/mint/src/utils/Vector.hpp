#ifndef MINT_UTILS_VECTOR_HXX_
#define MINT_UTILS_VECTOR_HXX_

#include "slic/slic.hpp"            // for slic macros

#include <cstring>                  // for std::memcpy
#include <cstdlib>                  // for std::malloc, std::realloc, st::free
#include <cmath>                    // for std::ceil

namespace axom {
namespace mint {

/*!
 * \class Vector
 * \brief Provides an array with dynamic reallocation and insertion.
 * \tparam T the type of the values to hold.
 * \tparam IndexType the integral type to use for indices.
 */
template< typename T, typename IndexType >
class Vector {
public:
  
  /*!
   * \brief Default constructor.
   */
  Vector():
    m_size(0),
    m_capacity(0),
    m_resize_ratio(2.0),
    m_chunk_size(1),
    m_data( AXOM_NULLPTR )
  {}

  /*!
   * \brief Custom constructor, creates a vector with the given capacity and
   *  resize ratio.
   * \param [in] capacity the initial capacity of the vector to allocate.
   * \param [in] ratio the scale factor by which to resize the vector when 
   *  size exceeds the capacity. 
   */
  Vector( IndexType capacity, double ratio=2.0 ):
    m_size(0),
    m_capacity( capacity ),
    m_resize_ratio( ratio ),
    m_chunk_size(1),
    m_data( AXOM_NULLPTR )
  {
    SLIC_ERROR_IF( m_capacity < 0, "Capacity cannot be negative." );
    if ( m_capacity > 0 ) {
      m_capacity = calculateCapacity( m_capacity );
      m_data = static_cast< T* >( std::malloc( m_capacity * sizeof(T) ) );
      SLIC_ERROR_IF( m_data == AXOM_NULLPTR, "Vector allocation failed." );
    } 
  }

  /*!
   * \brief Custom constructor, creates a vector with the given capacity,
   *  chunk size, and resize ratio.
   * \param [in] capacity the initial capacity of the vector to allocate.
   * \param [in] chunk_size when allocation is done it is done in chunks of this
   *  size.
   * \param [in] ratio the scale factor by which to resize the vector when 
   *  size exceeds the capacity. 
   */
  Vector( IndexType capacity, uint chunk_size, double ratio=2.0 ):
    m_size(0),
    m_capacity( capacity ),
    m_resize_ratio( ratio ),
    m_chunk_size( chunk_size ),
    m_data( AXOM_NULLPTR )
  {
    SLIC_ERROR_IF( m_capacity < 0, "Capacity cannot be negative." );
    if ( m_capacity > 0 ) {
      m_capacity = calculateCapacity( m_capacity );
      m_data = static_cast< T* >( std::malloc( m_capacity * sizeof(T) ) );
      SLIC_ERROR_IF( m_data == AXOM_NULLPTR, "Vector allocation failed." );
    } 
  }
  
  /*!
   * \brief Custom destructor, free's the data buffer.
   */
  ~Vector() {
    if ( m_data != AXOM_NULLPTR ) {
      std::free( m_data );
      m_data = AXOM_NULLPTR;
    }
  }

  /*!
   * \brief Accessor, returns a reference to the value at position pos.
   * \return a reference to the value at position pos.
   */
  inline T& operator[]( IndexType pos ) const
  { return m_data[ pos ]; }

  /*!
   * \brief Returns true iff the vector stores no elements.
   * \return true iff the vector stores no elements.
   * \note If the vector is empty the capacity can still be greater than zero.
   */
  inline bool empty() const 
  { return m_size == 0; }


  /*!
   * \brief Append a value to the end of the array.
   * \param [in] value the value to append.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void add( const T& value ) {
    SLIC_ASSERT( m_size <= m_capacity );

    if ( m_size + 1 > m_capacity ) {
      dynamicResize( m_size + 1 );
    }

    m_data[ m_size++ ] = value;
  }

  /*!
   * \brief Append values to the end of the array.
   * \param [in] values the values to append.
   * \param [in] n the number of values to append.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void add( const T* values, IndexType n ) {
    SLIC_ASSERT( m_size <= m_capacity );

    if ( m_size + n > m_capacity ) {
      dynamicResize( m_size + n );
    }
    
    std::memcpy( m_data + m_size, values, n * sizeof(T) );
    m_size += n;
  }

  /*!
   * \brief Write values into existing spots in the array.
   * \param [in] values the new values to be written.
   * \param [in] n the number of values to write.
   * \param [in] pos the position at which to begin writing.
   * \pre pos + n <= m_size.
   */
  inline void set( const T* values, IndexType n, IndexType pos ) const {
    SLIC_ASSERT( pos + n <= m_size );
    std::memcpy( m_data + pos, values, n * sizeof(T) );
  }

  /*!
   * \brief Make space for a subsequent insertion into the array.
   * \param [in] n the number of values to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline T* reserveForInsert( IndexType n, IndexType pos ) {
    SLIC_ASSERT( m_size <= m_capacity );

    if ( m_size + n > m_capacity ) {
      dynamicResize( m_size + n );
    }

    std::memcpy( m_data + pos + n, m_data + pos, ( m_size - pos ) * sizeof(T) );
    return m_data + pos;
  }

  /*!
   * \brief Insert values into the array at the given position.
   * \param [in] values the values to insert.
   * \param [in] n the number of values to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T* values, IndexType n, IndexType pos ) {
    SLIC_ASSERT( m_size <= m_capacity );
    reserveForInsert( n, pos );
    std::memcpy( m_data + pos, values, n * sizeof(T) );
    m_size += n;
  }

  /*!
   * \brief Return a pointer to the array of data.
   * \return a pointer to the array of data.
   */
  inline T* getData() const
  { return m_data; }

  /*!
   * \brief Return the amount of space allocated for the data array.
   * \return the amount of space allocated for the data array.
   */
  inline IndexType getCapacity() const
  { return m_capacity; }

  /*!
   * \brief Set the ammount of space allocated for the data array.
   * \param [in] capacity the new amount of space to allocate.
   * \note the capacity of the array is rounded up to the nearest multiple of
   *  m_chunk_size.
   */
  inline void setCapacity( IndexType capacity ) {
    SLIC_ASSERT( capacity >= 0 );

    m_capacity = calculateCapacity( capacity );
    m_data = static_cast< T* >( std::realloc( m_data, m_capacity * sizeof(T) ) );
    SLIC_ERROR_IF( m_data == AXOM_NULLPTR && capacity > 0, 
                   "Vector reallocation failed." );

    if ( m_size > m_capacity ) {
      m_size = m_capacity;
    }
  }

  /*!
   * \brief Return the number of values stored in the data array.
   * \return the number of values stored in the data array.
   */
  inline IndexType getSize() const
  { return m_size; }


  /*!
   * \brief Set the number of values stored in the data array.
   * \param [in] size the new size of the data array.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void setSize( IndexType size ) {
    SLIC_ASSERT( size >= 0 );
    if ( size > m_capacity ) {
      setCapacity( size );
    }

    m_size = size;
  }

 /*!
  * \brief Get the ratio by which the capacity increases upon dynamic resize. 
  * \return The ration by which the capacity increases upon dynamic resize.
  */
  inline double getResizeRatio() const
  { return m_resize_ratio; }
 
 /*!
  * \brief Set the ratio by which the capacity increases upon dynamic resize. 
  * \param [in] ratio the scale factor by which to resize the vector when 
  *  size exceeds the capacity.
  */
  inline void setResizeRatio( double ratio ) 
  { m_resize_ratio = ratio; }

  /*!
  * \brief Get the chunk size of all allocations. 
  * \return the chunk size of all allocations.
  */
  inline uint getChunkSize() const
  { return m_chunk_size; }

  /*!
  * \brief Set the chunk size of all allocations. 
  * \param [in] chunk_size when allocation is done it is done in chunks of this
  *  size.
  */
  inline void setChunkSize( uint chunk_size )
  { m_chunk_size = chunk_size; }
  

private:

  /*!
  * \brief Returns the smallest multiple of m_chunk_size greater than or equal
  *  to the given capacity. 
  * \param [in] capacity the capacity to round up.
  * \return The capacity rounded up to the nearest multiple of m_chunk_size.
  */
  inline IndexType calculateCapacity( IndexType capacity ) {
    return std::ceil( double( capacity ) / m_chunk_size ) * m_chunk_size;
  }

  /*!
  * \brief Reallocates the data array when the size exceeds the capacity. 
  * \param [in] newSize the size of the data array which exceeds the current
  *  capacity.
  */
  inline void dynamicResize( IndexType newSize ) {
    SLIC_ERROR_IF( m_resize_ratio < 1.0, "Resize ratio of " << m_resize_ratio <<
                                         " doesn't support dynamic resizing");

    m_capacity = calculateCapacity( newSize * m_resize_ratio );
    m_data = static_cast< T* >( std::realloc( m_data, m_capacity * sizeof(T) ) );
    SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0, 
                   "Vector reallocation failed." );
  }

  IndexType m_size;
  IndexType m_capacity;
  double m_resize_ratio;
  uint m_chunk_size;
  T* m_data;
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_UTILS_VECTOR_HXX_ */
