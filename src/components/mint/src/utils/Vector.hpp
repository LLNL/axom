#ifndef MINT_UTILS_VECTOR_HXX_
#define MINT_UTILS_VECTOR_HXX_

#include "slic/slic.hpp"                // for slic macros
#include "mint/DataTypes.hpp"
#include "axom_utils/Utilities.hpp"

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

#include <cstring>                      // for std::memcpy
#include <cmath>                        // for std::ceil

namespace axom
{
namespace mint
{
/*!
 * \class Vector
 * \brief Provides an array with dynamic reallocation and insertion.
 * \tparam T the type of the values to hold.
 */
template< typename T >
class Vector
{
private:
  static constexpr double DEFAULT_RESIZE_RATIO=2.0;

public:

  /*!
   * \brief Default constructor.
   */
  Vector() = delete;

  /*!
   * \brief Constructs a Vector instance with the give number of tuples.
   *
   * \param [in] num_tuples the number of tuples the vector holds.
   * \param [in] num_components the number of components per tuple.
   */
  Vector( localIndex num_tuples, localIndex num_components=1 ) :
#ifdef MINT_USE_SIDRE
   m_view( AXOM_NULLPTR ),
#endif
   m_data( AXOM_NULLPTR ),
   m_num_tuples( num_tuples ),
   m_capacity( DEFAULT_RESIZE_RATIO*num_tuples ),
   m_num_components( num_components ),
   m_resize_ratio( DEFAULT_RESIZE_RATIO )
  {

  }

  /*!
   * \brief Constructor for use without sidre.
   * \param [in] num_components the number of values per tuple.
   * \param [in] capacity the number of tuples to allocate space for.
   * \param [in] num_tuples the number of tuples accounted for in the vector.
   * \param [in] ratio the ratio by which to resize the vector when the
   *  size exceeds the capacity. A ratio less than one prohibits resizing.
   */
  Vector( localIndex num_components,
          localIndex capacity,
          localIndex num_tuples,
          double ratio ) :
#ifdef MINT_USE_SIDRE
    m_view( AXOM_NULLPTR ),
#endif
    m_data( AXOM_NULLPTR ),
    m_num_tuples( num_tuples ),
    m_capacity( capacity ),
    m_num_components( num_components ),
    m_resize_ratio( ratio )
  {
    SLIC_ERROR_IF( m_capacity < 0,
                   "Tuple capacity (" << m_capacity << ") " <<
                   "cannot be negative." );
    SLIC_ERROR_IF( m_num_tuples < 0,
                   "Number of tuples (" << m_num_tuples << ") " <<
                   "cannot be negative." );
    SLIC_ERROR_IF( m_num_components <= 0,
                   "Components per tuple (" << m_num_components << ") " <<
                   "must be greater than zero." );
    SLIC_ERROR_IF( m_num_tuples > m_capacity,
                   "Number of tuples (" << m_num_tuples << ") " <<
                   "cannot be greater than the tuple capacity " <<
                   "(" << m_capacity << ")." );

//    TODO: missing function ?
//    setCapacity( capacity );
  }

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Constructor for use with a sidre::View that already has data.
   * \param [in] view the sidre::View that holds this Vector's data.
   * \param [in] num_tuples the number of tuples accounted for in the vector.
   * \param [in] ratio the ratio by which to resize the vector when the
   *  size exceeds the capacity. A ratio less than one prohibits resizing.
   * \pre view != AXOM_NULLPTR.
   */
  Vector( sidre::View * view, const localIndex & capacity,
          const localIndex & num_tuples, const double & ratio ) :
    m_view( view ),
    m_data( AXOM_NULLPTR )
    m_num_tuples( num_tuples ),
    m_tuple_capacity( capacity ),
    m_tuple_components( getViewShape( view, 1 ) ),
    m_resize_ratio( ratio )
  {
    localIndex view_capacity = getViewShape( view, 0 );
    SLIC_ERROR_IF( view_capacity != m_tuple_capacity,
                   "View capacity (" << view_capacity << ") " <<
                   "does not equal the given capacity " <<
                   "(" << m_tuple_capacity << ")." );
    SLIC_ERROR_IF( m_tuple_capacity < 0,
                   "Tuple capacity (" << m_tuple_capacity << ") " <<
                   "cannot be negative." );
    SLIC_ERROR_IF( m_num_tuples < 0,
                   "Number of tuples (" << m_num_tuples << ") " <<
                   "cannot be negative." );
    SLIC_ERROR_IF( m_tuple_components <= 0,
                   "Components per tuple (" << m_tuple_components << ") " <<
                   "must be greater than zero." );
    SLIC_ERROR_IF( m_num_tuples > m_tuple_capacity,
                   "Number of tuples (" << m_num_tuples << ") " <<
                   "cannot be greater than the tuple capacity " <<
                   "(" << m_tuple_capacity << ")." );

    sidre::TypeID view_type = view->getTypeID();
    sidre::TypeID T_type = sidre::detail::SidreTT< T >::id;
    SLIC_ERROR_IF( view_type != T_type, "View data type (" << view_type << ")"
                      << "differs from this vector type (" << T_type << ")." );

    m_data = view->getVoidPtr();
  }

  /*!
   * \brief Constructor for use with an empty sidre::View which will hold data.
   * \param [in] view the sidre::View that will holds this Vector's data.
   * \param [in] num_components the number of values per tuple.
   * \param [in] capacity the number of tuples to allocate space for.
   * \param [in] num_tuples the number of tuples accounted for in the vector.
   * \param [in] ratio the ratio by which to resize the vector when the
   *  size exceeds the capacity. A ratio less than one prohibits resizing.
   * \pre view != AXOM_NULLPTR.
   */
  Vector( sidre::View * view, localIndex num_components,
          const localIndex & capacity, const localIndex & num_tuples,
          const double & ratio ) :
    m_view( view ),
    m_data( AXOM_NULLPTR )
    m_num_tuples( num_tuples ),
    m_tuple_capacity( capacity ),
    m_tuple_components( num_components ),
    m_resize_ratio( ratio ),
  {
    SLIC_ERROR_IF( m_tuple_capacity < 0,
                   "Tuple capacity (" << m_tuple_capacity << ") " <<
                   "cannot be negative." );
    SLIC_ERROR_IF( m_num_tuples < 0,
                   "Number of tuples (" << m_num_tuples << ") " <<
                   "cannot be negative." );
    SLIC_ERROR_IF( m_tuple_components <= 0,
                   "Components per tuple (" << m_tuple_components << ") " <<
                   "must be greater than zero." );
    SLIC_ERROR_IF( m_num_tuples > m_tuple_capacity,
                   "Number of tuples (" << m_num_tuples << ") " <<
                   "cannot be greater than the tuple capacity " <<
                   "(" << m_tuple_capacity << ")." );
    SLIC_ERROR_IF( !m_view->isEmpty, "View must be empty." );

    sidre::TypeID T_type = sidre::detail::SidreTT< T >::id;
    sidre::SidreLength dims[ 2 ];
    dims[0] = m_tuple_capacity;
    dims[1] = m_tuple_components;
    view->describe( T_type, 2, dims );
    view->allocate();
    m_data = view->getVoidPtr();
  }
#endif

  /*!
   * \brief Accessor, returns a reference to the value at position pos.
   * \return a reference to the value at position pos.
   */
  inline T & operator[]( localIndex pos ) { return m_data[ pos ]; }

  /*!
   * \brief Constant accessor, returns a constant reference to the value at
   *  position pos.
   * \return a reference to the value at position pos.
   */
  inline const T & operator[]( localIndex pos ) const { return m_data[ pos ]; }

  /*!
   * \brief Returns true iff the vector stores no elements.
   * \return true iff the vector stores no elements.
   * \note If the vector is empty the capacity can still be greater than zero.
   */
  inline bool empty() const { return m_num_tuples == 0; }

  /*!
   * \brief Append a value to the end of the array.
   * \param [in] value the value to append.
   * \pre m_tuple_components == 1.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void add( const T& value )
  {
    SLIC_ASSERT_MSG( m_num_components != 1,
                    "Number of components must be 1." );

    localIndex new_size = m_num_tuples + 1;
    if ( new_size > m_capacity )
    { dynamicRealloc( new_size ); }

    m_data[ m_num_tuples ] = value;
  }

  /*!
   * \brief Append values to the end of the array.
   * \param [in] values the values to append.
   * \param [in] n the number of values to append.
   * \pre n % m_num_tuples == 0.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void add( const T * values, localIndex n )
  {
    SLIC_ASSERT_MSG( n % m_num_components != 0, "n (" << n << ") must be a "
                      << "multiple of the number of components of this vector ("
                      << m_num_components << ")." );

    localIndex new_size = m_num_tuples + ( n / m_num_components );
    if ( new_size > m_capacity )
    { dynamicRealloc( new_size ); }

    T * cur_end = m_data + m_num_tuples * m_num_components;
    std::memcpy( m_data , values, n * sizeof(T) );
  }

  /*!
   * \brief Write values into existing spots in the array.
   * \param [in] values the new values to be written.
   * \param [in] n the number of values to write.
   * \param [in] pos the position at which to begin writing.
   * \pre pos + n <= m_num_tuples.
   */
  inline void set( const T * values, localIndex n, localIndex pos )
  {
    SLIC_ASSERT( pos >= 0 );
    SLIC_ASSERT( pos + n <= m_num_components * m_num_tuples );
    std::memcpy( m_data + pos, values, n * sizeof(T) );
  }

  /*!
   * \brief Make space for a subsequent insertion into the array.
   * \param [in] n the number of values to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \pre n % m_num_tuples == 0.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline T * reserveForInsert( localIndex n, localIndex pos )
  {
    SLIC_ASSERT( pos <= m_num_tuples * m_num_components );
    SLIC_ASSERT_MSG( n % m_num_components != 0, "n (" << n << ") must be a "
                      << "multiple of the number of components of this vector ("
                      << m_num_components << ")." );

    localIndex new_size = m_num_tuples + ( n / m_num_components );
    if ( new_size > m_capacity )
    {
      dynamicRealloc( new_size );
    }

    const T * const insert_pos = m_data + pos;
    T * cur_pos = m_data + m_num_tuples - 1;
    for ( ; cur_pos >= insert_pos ; --cur_pos )
    { *(cur_pos + n) = *cur_pos; }

    return m_data + pos;
  }

  /*!
   * \brief Insert values into the array at the given position.
   * \param [in] values the values to insert.
   * \param [in] n the number of values to insert.
   * \param [in] pos the position at which to begin the insertion.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void insert( const T * values, localIndex n, localIndex pos )
  {
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
    reserveForInsert( 1, pos );
    m_data[ pos ] = value;
  }

  /*!
   * \brief Return a pointer to the array of data.
   * \return a pointer to the array of data.
   */
  inline T * getData() { return m_data; }

  /*!
   * \brief Return a constant pointer to the array of data.
   * \return a pointer to the array of data.
   */
  inline const T * getData() const { return m_data; }

  /*!
   * \brief Return the number of tuples allocated for the data array.
   * \return the amount of space allocated for the data array.
   */
  inline localIndex getCapacity() const { return m_capacity; }

  /*!
   * \brief Set the number of tuples allocated for the data array.
   * \param [in] capacity the new number of tuples to allocate.
   */
  inline void updateTupleCapacity()
  {
    SLIC_ASSERT( m_capacity >= 0 );

#ifdef MINT_USE_SIDRE
    if ( m_view != AXOM_NULLPTR )
    { return allocate_view_data(); }
#endif

    m_data = utilities::realloc( m_data, m_capacity );
    SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                   "Vector reallocation failed." );
   }

  /*!
   * \brief Return the number of tuples stored in the data array.
   * \return the number of tuples stored in the data array.
   */
  inline localIndex getNumTuples() const { return m_num_tuples; }

  /*!
   * \brief Update the number of tuples stored in the data array.
   * \note Reallocation is done if the new size will exceed the capacity.
   */
  inline void updateNumTuples()
  {
    SLIC_ASSERT( m_num_tuples >= 0 );

    if ( m_num_tuples > m_capacity )
    {
//      TODO: missing function?
//      setCapacity( num_tuples );
    }

  }

  /*!
   * \brief Get the ratio by which the capacity increases upon dynamic resize.
   * \return The ration by which the capacity increases upon dynamic resize.
   */
  inline double getResizeRatio() const { return m_resize_ratio; }

  /*!
   * \brief Get the chunk size of all allocations.
   * \return the chunk size of all allocations.
   */
  inline localIndex getNumComponents() const { return m_num_components; }

private:

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Given a non-empty sidre::View of dimension 2, returns the length
   *  of the given dimension.
   * \param [in] view the sidre::View to examine.
   * \param [in] dim the dimension (0 or 1) to return the length of.
   * \return The length of dimension dim of view.
   */
  inline localIndex getViewShape( const sidre::View * view, bool dim ) const
  {
    SLIC_ERROR_IF( view == AXOM_NULLPTR, "view cannot be a null pointer." );
    SLIC_ERROR_IF( view->isEmpty(), "view cannot be empty." );
    SLIC_ERROR_IF( view->getNumDimensions() != 2,
                                                 "view must have dimension 2.");

    sidre::SidreLength dims[ 2 ];
    view->getShape( 2, dims );
    return dims[ dim ];
  }

  /*!
   * \brief Allocates space within the Vector's sidre::View.
   */
  inline void allocate_view_data()
  {
    view->reallocate( m_tuple_capacity * m_tuple_components );

    sidre::TypeID T_type = sidre::detail::SidreTT< T >::id;
    sidre::SidreLength dims[ 2 ];
    dims[0] = m_tuple_capacity;
    dims[1] = m_tuple_components;
    view->describe( T_type, 2, dims );
    view->apply();
    m_data = view->getVoidPtr();

    SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_tuple_capacity > 0,
                   "Vector reallocation failed." );
  }
#endif

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   * \param [in] new_n_tuples the number of tuples which exceeds the current
   *  capacity.
   */
  inline void dynamicRealloc( localIndex new_n_tuples )
  {
    SLIC_ERROR_IF( m_resize_ratio < 1.0, "Resize ratio of " << m_resize_ratio <<
                   " doesn't support dynamic resizing");
    m_capacity = new_n_tuples * m_resize_ratio + 0.5;

#ifdef MINT_USE_SIDRE
    if ( m_view != AXOM_NULLPTR )
    { return allocate_view_data(); }
#endif

    m_data = utilities::realloc( m_data, m_capacity );

    SLIC_ERROR_IF( m_data == AXOM_NULLPTR && m_capacity > 0,
                   "Vector reallocation failed." );
  }

#ifdef MINT_USE_SIDRE
  sidre::View * m_view
#endif

  T* m_data;
  localIndex m_num_tuples;
  localIndex m_capacity;
  localIndex m_num_components;
  double m_resize_ratio;

};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_UTILS_VECTOR_HXX_ */
