#ifndef MINT_UTILS_VECTOR_HXX_
#define MINT_UTILS_VECTOR_HXX_

#include "slic/slic.hpp"                // for slic macros
#include "mint/DataTypes.hpp"
#include "axom_utils/Utilities.hpp"


#include <cstring>                      // for std::memcpy
#include <cmath>                        // for std::ceil

namespace axom {
namespace mint {

/*!
 * \class Vector
 * \brief Provides an array with dynamic reallocation and insertion.
 * \tparam T the type of the values to hold.
 */
template< typename T >
class Vector 
{
public:
  
  /*!
   * \brief Default constructor.
   */
  Vector():
    m_size(0),
    m_capacity(0),
    m_resize_ratio(2.0),
    m_components(1),
    m_data( AXOM_NULLPTR )
  {}

  /*!
   * \brief Custom constructor, creates a vector with the given capacity and
   *  resize ratio.
   * \param [in] capacity the initial capacity of the vector to allocate.
   * \param [in] ratio the scale factor by which to resize the vector when 
   *  size exceeds the capacity. 
   */
  Vector( localIndex capacity, double ratio=2.0 ):
    m_size(0),
    m_resize_ratio( ratio ),
    m_components(1),
    m_data( AXOM_NULLPTR )
  {
    setCapacity( capacity );
  }

  /*!
   * \brief Custom constructor, creates a vector with the given capacity,
   *  chunk size, and resize ratio.
   * \param [in] capacity the initial capacity of the vector to allocate.
   * \param [in] num_components the number of components per value. When 
   *  allocation is done the resulting capacity is a multiples of the number
   *  of components.
   * \param [in] ratio the scale factor by which to resize the vector when 
   *  size exceeds the capacity. 
   */
  Vector( localIndex capacity, int num_components, double ratio=2.0 ):
    m_size(0),
    m_resize_ratio( ratio ),
    m_components( num_components ),
    m_data( AXOM_NULLPTR )
  {
    setCapacity( capacity );
  }

  /*!
   * \breif Copy constructor.
   * \param [in] rhs the Vector to copy.
   */
  Vector( const Vector & rhs ):
    m_size(rhs.m_size),
    m_resize_ratio(rhs.m_resize_ratio),
    m_components(rhs.m_components)
  {
    setCapacity( rhs.m_capacity );
    localIndex byte_size = m_size * sizeof(T); 
    std::memcpy(m_data, rhs.m_data, byte_size);
  }

  /*!
   * \breif Move constructor.
   * \param [in] rhs the Vector to move.
   */
  Vector( const Vector && rhs ):
    m_size(rhs.m_size),
    m_resize_ratio(rhs.m_resize_ratio),
    m_components(rhs.m_components),
    m_data(rhs.m_data)
  {
    setCapacity( rhs.m_capacity );

    if (&rhs != this) {
      rhs.m_data = AXOM_NULLPTR;    
    }
  }
  
  /*!
   * \brief Custom destructor, free's the data buffer.
   */
  ~Vector() 
  {
    utilities::free( m_data );
    m_data = AXOM_NULLPTR;
  }

  /*!
   * \breif Copy assignment operator.
   * \param [in] rhs the Vector to copy.
   */
  Vector< T > & operator=( const Vector & rhs )
  { 
    m_size = rhs.m_size;
    setCapacity( rhs.m_capacity );
    m_resize_ratio = rhs.m_resize_ratio;
    m_components = rhs.m_components;
    localIndex byte_size = m_size * sizeof(T); 
    std::memcpy(m_data, rhs.m_data, byte_size);
    return *this;
  }

  /*!
   * \breif Move assignment operator.
   * \param [in] rhs the Vector to move.
   */
  Vector< T > & operator=( const Vector && rhs )
  { 
    m_size = rhs.m_size;
    setCapacity( rhs.m_capacity );
    m_resize_ratio = rhs.m_resize_ratio;
    m_components = rhs.m_components;
    m_data = rhs.m_data;

    if (&rhs != this) {
      rhs.m_data = AXOM_NULLPTR;
    }
  }

  /*!
   * \brief Accessor, returns a reference to the value at position pos.
   * \return a reference to the value at position pos.
   */
  inline T& operator[]( localIndex pos ) const
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
  inline void add( const T& value ) 
  {
    SLIC_ASSERT( m_size <= m_capacity );

    if ( m_size + 1 > m_capacity ) {
      dynamicRealloc( m_size + 1 );
    }

    m_data[ m_size++ ] = value;
  }

  /*!
   * \brief Append values to the end of the array.
   * \param [in] values the values to append.
   * \param [in] n the number of values to append.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void add( const T* values, localIndex n ) 
  {
    SLIC_ASSERT( m_size <= m_capacity );

    if ( m_size + n > m_capacity ) {
      dynamicRealloc( m_size + n );
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
  inline void set( const T* values, localIndex n, localIndex pos ) const 
  {
    SLIC_ASSERT( m_size <= m_capacity );
    SLIC_ASSERT( pos >= 0 && pos + n <= m_size );
    std::memcpy( m_data + pos, values, n * sizeof(T) );
  }

  /*!
   * \brief Make space for a subsequent insertion into the array.
   * \param [in] n the number of values to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline T* reserveForInsert( localIndex n, localIndex pos ) 
  {
    SLIC_ASSERT( m_size <= m_capacity );
    SLIC_ASSERT( pos <= m_size );

    if ( m_size + n > m_capacity ) {
      dynamicRealloc( m_size + n );
    }

    const T* const insert_pos = m_data + pos;
    T* cur_pos = m_data + m_size - 1;
    for (; cur_pos >= insert_pos; --cur_pos ) {
      *(cur_pos + n) = *cur_pos;
    }

    m_size += n;
    return m_data + pos;
  }

  /*!
   * \brief Insert values into the array at the given position.
   * \param [in] values the values to insert.
   * \param [in] n the number of values to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T* values, localIndex n, localIndex pos ) 
  {
    SLIC_ASSERT( m_size <= m_capacity );

    reserveForInsert( n, pos );
    std::memcpy( m_data + pos, values, n * sizeof(T) );
  }

  /*!
   * \brief Insert a value into the array at the given position.
   * \param [in] value the value to insert.
   * \param [in] pos the position at which to insert.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T& value, localIndex pos )
  {
    SLIC_ASSERT( m_size <= m_capacity );

    reserveForInsert( 1, pos );
    m_data[ pos ] = value;
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
  inline localIndex getCapacity() const
  { return m_capacity; }

  /*!
   * \brief Set the ammount of space allocated for the data array.
   * \param [in] capacity the new amount of space to allocate.
   * \note the capacity of the array is rounded up to the nearest multiple of
   *  m_components.
   */
  inline void setCapacity( localIndex capacity ) 
  {
    SLIC_ASSERT( capacity >= 0 );

    m_capacity = std::ceil( 1.0 * capacity / m_components ) * m_components;
    m_data = utilities::realloc( m_data, m_capacity );
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
  inline localIndex getSize() const
  { return m_size; }


  /*!
   * \brief Set the number of values stored in the data array.
   * \param [in] size the new size of the data array.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void setSize( localIndex size ) 
  {
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
  inline int getNumComponents() const
  { return m_components; }

  /*!
  * \brief Get the number of tuples. 
  * \return the number of tuples.
  */
  inline localIndex getNumTuples() const
  {
    SLIC_WARNING_IF( m_size % m_components != 0, "Size of Vector is not a " <<
                     "multiple of the number of components.");
    return m_size / m_components;
  }
  

private:

  /*!
  * \brief Reallocates the data array when the size exceeds the capacity. 
  * \param [in] newSize the size of the data array which exceeds the current
  *  capacity.
  */
  inline void dynamicRealloc( localIndex newSize ) 
  {
    SLIC_ERROR_IF( m_resize_ratio < 1.0, "Resize ratio of " << m_resize_ratio <<
                                         " doesn't support dynamic resizing");
    m_capacity = std::ceil( newSize * m_resize_ratio / m_components ) * m_components;

    m_data = utilities::realloc( m_data, m_capacity );

#ifdef USE_SIDRE
  if ( m_view != AXOM_NULLPTR ) {
    m_data = m_view->reallocate()->getPointer();
  } else {
    m_data = utilities::realloc( m_data, m_capacity );
  }
#endif    

    SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0, 
                   "Vector reallocation failed." );
  }

  localIndex m_size;
  localIndex m_capacity;
  double m_resize_ratio;
  int m_components;
  T* m_data;

#ifdef USE_SIDRE
  sidre::View * m_view
#endif
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_UTILS_VECTOR_HXX_ */
