/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "mint/ConnectivityArray.hpp"
#include "mint/CellTypes.hpp"
#include "mint/config.hpp"
#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"      /* for UnitTestLogger */

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#endif

#include <algorithm>

const char IGNORE_OUTPUT[] = ".*";

namespace axom
{
namespace mint
{
namespace internal
{


/*!
 * \brief Return the new ID capacity for a ConnectivityArray given an 
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the ammount the number of IDs will increase by
 *
 * \note this template specialization is for the ConnectivityArray with fixed
 *  stride.
 */
template< IndexType STRIDE, CellTypes TYPE >
IndexType calc_ID_capacity( const ConnectivityArray< STRIDE, TYPE >& connec, 
                            IndexType increase )
{
  IndexType new_n_IDs = connec.getNumberOfIDs() + increase;
  if ( new_n_IDs > connec.getIDCapacity() )
  {
    return new_n_IDs * connec.getResizeRatio() + 0.5;
  }

  return connec.getIDCapacity();
}

/*!
 * \brief Return the new ID capacity for a ConnectivityArray given an 
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the ammount the number of IDs will increase by
 *
 * \note this template specialization is for the ConnectivityArray with
 *  indirection but fixed type. The plus/minus one are to account for the
 *  fact that the offset array which has capacity getIDCapacity() + 1.
 */
template< CellTypes TYPE >
IndexType calc_ID_capacity( const ConnectivityArray< NEED_INDIRECTION, TYPE >& connec, 
                            IndexType increase )
{
  IndexType new_n_IDs = connec.getNumberOfIDs() + increase;
  if ( new_n_IDs > connec.getIDCapacity() )
  {
    return ((new_n_IDs + 1) * connec.getResizeRatio() + 0.5) - 1;
  }

  return connec.getIDCapacity();
}

/*!
 * \brief Return the new ID capacity for a ConnectivityArray given an 
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the ammount the number of IDs will increase by
 *
 * \note this template specialization is for the ConnectivityArray with
 *  indirection and variable type. No plus/minus one is needed since the
 *  types Array is used to calculate ID capacity which has capacity
 *  getIDCapacity().
 */
IndexType calc_ID_capacity( const ConnectivityArray< NEED_INDIRECTION, CellTypes::MIXED >& connec, 
                            IndexType increase )
{
  IndexType new_n_IDs = connec.getNumberOfIDs() + increase;
  if ( new_n_IDs > connec.getIDCapacity() )
  {
    return new_n_IDs * connec.getResizeRatio() + 0.5;
  }

  return connec.getIDCapacity();
}

/*!
 * \brief Return the new value capacity for a ConnectivityArray given an 
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the ammount the number of values will increase by
 *
 * \note this template specialization is for the ConnectivityArray with fixed
 *  stride.
 */
template< IndexType STRIDE, CellTypes TYPE >
IndexType calc_value_capacity( const ConnectivityArray< STRIDE, TYPE >& connec, 
                              IndexType increase )
{ return calc_ID_capacity( connec, increase / STRIDE ) * STRIDE; }


/*!
 * \brief Return the new value capacity for a ConnectivityArray given an 
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the ammount the number of values will increase by
 *
 * \note this template specialization is for the ConnectivityArray with
 *  indirection.
 */
template< CellTypes TYPE >
IndexType calc_value_capacity( const ConnectivityArray< NEED_INDIRECTION, TYPE >& connec, 
                              IndexType increase )
{
  IndexType new_n_values = connec.getNumberOfValues() + increase;
  if ( new_n_values > connec.getValueCapacity() )
  {
    return new_n_values * connec.getResizeRatio() + 0.5;
  }

  return connec.getValueCapacity();
}

/*!
 * \brief Get the type associated with the given ID to be inserted into a
 *  ConnectivityArray.
 *
 * \param [in] ID the ID in question.
 * \param [in] types the array of types.
 *
 * \tparam TYPE the type of the ConnectivityArray to be inseted into.
 */
template< CellTypes TYPE >
CellTypes get_type( IndexType ID, const CellTypes* types )
{
  if ( types == AXOM_NULLPTR )
  {
    return TYPE;
  }
  else
  {
    return types[ ID ];
  }
}

/*!
 * \brief Get a pointer to the types starting at the given ID.
 *
 * \param [in] ID the ID in question.
 * \param [in] types the array of types.
 *
 * \tparam TYPE the type of the ConnectivityArray to be inseted into.
 */
const CellTypes* get_type_ptr( IndexType ID, const CellTypes* types )
{
  if ( types == AXOM_NULLPTR )
  {
    return AXOM_NULLPTR;
  }
  else
  {
    return types +  ID;
  }
}

/*!
 * \brief Return the number of values associated with ID.
 *
 * \param [in] ID the ID in question.
 * \param [in] offsets the offsets array.
 *
 * \tparam STRIDE the stride of the ConnectivityArray to be inseted into.
 */
template< IndexType STRIDE  >
IndexType get_n_values( IndexType ID, const IndexType* offsets )
{
  if ( offsets == AXOM_NULLPTR )
  {
    return STRIDE;
  }
  else
  {
    return offsets[ ID + 1 ] - offsets[ ID ];
  }
}

/*!
 * \brief Return the number of values between the given IDs.
 *
 * \param [in] startID the start of the range.
 * \param [in] endID the end of the range.
 *
 * \tparam STRIDE the stride of the ConnectivityArray to be inseted into.
 */
template< IndexType STRIDE  >
IndexType get_n_values( IndexType startID, IndexType endID, 
                        const IndexType* offsets )
{
  if ( offsets == AXOM_NULLPTR )
  {
    return STRIDE * (endID - startID);
  }
  else
  {
    return offsets[ endID ] - offsets[ startID ];
  }
}

/*!
 * \brief Return a pointer to the values array starting with the given ID.
 *
 * \param [in] ID the ID in question.
 * \param [in] values the values array.
 * \param [in] offsets the offsets array.
 *
 * \tparam STRIDE the stride of the ConnectivityArray to be inseted into.
 */
template< IndexType STRIDE >
const IndexType* get_value_ptr( IndexType ID, const IndexType* values, 
                                const IndexType* offsets  )
{
  if ( offsets == AXOM_NULLPTR )
  {
    return values + ID * STRIDE;
  }
  else
  {
    return values + offsets[ ID ];
  }
}

/*!
 * \brief Return a pointer to the offsets array starting with the given ID.
 *
 * \param [in] ID the ID in question.
 * \param [in] offsets the offsets array.
 */
const IndexType* get_offset_ptr( IndexType ID, const IndexType* offsets  )
{
  if ( offsets == AXOM_NULLPTR )
  {
    return AXOM_NULLPTR;
  }
  else
  {
    return offsets + ID;
  }
}

/*!
 * \brief Check that the given pointers point to the same values or that they
 *  are both null.
 *
 * \param [in] p1 the first pointer to check.
 * \param [in] p2 the second pointer to check.
 * \param [in] n the number of values to check.
 * \param [in] line the line from which this function was called.
 */
template < typename T >
void check_pointers( const T* p1, const T* p2, IndexType n, int line )
{
  if ( p1 == AXOM_NULLPTR && p2 == AXOM_NULLPTR )
  {
    return;
  }

  for ( IndexType i = 0; i < n; ++i )
  {
    EXPECT_EQ( p1[ i ], p2[ i ] ) << "i = " << i << ". Called from line " 
                                                               << line << ".\n";
  }
}

/*!
 * \brief Check that the given ConnectivityArrays are equal.
 *
 * \param [in] lhs the first ConnectivityArray to check.
 * \param [in] rhs the second ConnectivityArray to check.
 * \param [in] line the line from which this function was called.
 *
 * \note does not compare resize ratios.
 */
template< IndexType STRIDE, CellTypes TYPE >
void check_equality( const ConnectivityArray< STRIDE, TYPE >& lhs, 
                     const ConnectivityArray< STRIDE, TYPE >& rhs, int line )
{
  ASSERT_EQ( lhs.getNumberOfIDs(), rhs.getNumberOfIDs() );
  ASSERT_EQ( lhs.getNumberOfValues(), rhs.getNumberOfValues() );
  EXPECT_EQ( lhs.getIDCapacity(), rhs.getIDCapacity() );
  EXPECT_EQ( lhs.isExternal(), rhs.isExternal() );
  EXPECT_EQ( lhs.isInSidre(), rhs.isInSidre() );

  const IndexType n_IDs = lhs.getNumberOfIDs();
  for ( IndexType ID = 0; ID < n_IDs; ++ID )
  {
    const IndexType cur_n_values = lhs.getNumberOfValuesForID( ID );
    ASSERT_EQ( cur_n_values, rhs.getNumberOfValuesForID( ID ) );
    EXPECT_EQ( lhs.getIDType( ID ), rhs.getIDType( ID ) );
    check_pointers( lhs[ ID ], rhs[ ID ], cur_n_values, line );
  }
}

/*!
 * \brief Check that appending to the given ConnectivityArray functions as 
 *  expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to.
 * \param [in] n_IDs the number of IDs to append.
 * \param [in] values the array of values to append.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 *
 * \note offsets only need to be specified for ConnectivityArrays that need 
 *  indirection and types only needs to be specified for ConnectivityArrays
 *  that have variable ID types. First appends half of the IDs one at a time and
 *  then appends the other half all at once.
 */
template< IndexType STRIDE, CellTypes TYPE >
void append( ConnectivityArray< STRIDE, TYPE >& connec, IndexType n_IDs,
             const IndexType* values, const IndexType* offsets=AXOM_NULLPTR, 
             const CellTypes* types=AXOM_NULLPTR )
{
  const IndexType initial_n_IDs = connec.getNumberOfIDs();
  IndexType half_n_IDs = n_IDs / 2;
  for ( IndexType ID = 0; ID < half_n_IDs; ++ID )
  {
    const IndexType cur_n_values = get_n_values< STRIDE >( ID, offsets );
    const IndexType* cur_values = get_value_ptr< STRIDE >( ID, values, offsets );
    const CellTypes type = get_type< TYPE >( ID, types );

    connec.append( cur_values, cur_n_values, type );
    EXPECT_EQ( connec.getNumberOfIDs(), initial_n_IDs + ID + 1 );
    EXPECT_EQ( connec.getIDType( initial_n_IDs + ID ), type );
    EXPECT_EQ( connec.getNumberOfValuesForID( initial_n_IDs + ID ), cur_n_values );
    check_pointers( cur_values, connec[ initial_n_IDs + ID ], cur_n_values, 
                    __LINE__ );
  }

  const IndexType cur_ID = half_n_IDs;
  const IndexType* cur_values = get_value_ptr< STRIDE >( cur_ID, values, offsets );
  const IndexType* cur_offsets = get_offset_ptr( cur_ID, offsets );
  const CellTypes* cur_types = get_type_ptr( cur_ID, types );
  connec.appendM( cur_values, n_IDs - cur_ID, cur_offsets, cur_types );
  
  EXPECT_EQ( connec.getNumberOfIDs(), initial_n_IDs + n_IDs );
  for ( IndexType ID = 0; ID < n_IDs; ++ID )
  {
    const IndexType cur_n_values = get_n_values< STRIDE >( ID, offsets);
    const IndexType* cur_values = get_value_ptr< STRIDE >( ID, values, offsets );
    const CellTypes cur_type = get_type< TYPE >( ID, types );
    EXPECT_EQ( connec.getIDType( initial_n_IDs + ID ), cur_type );
    EXPECT_EQ( connec.getNumberOfValuesForID( initial_n_IDs + ID ), cur_n_values );
    check_pointers( cur_values, connec[ initial_n_IDs + ID ], cur_n_values, 
                    __LINE__ );
  }
}

/*!
 * \brief Check that setting the values of the given ConnectivityArray 
 *  functions as expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to.
 * \param [in] n_IDs the number of IDs to append.
 * \param [in] initial_values the array of values to append.
 * \param [in] values the array of values to set.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 *
 * \note offsets only need to be specified for ConnectivityArrays that need 
 *  indirection and types only needs to be specified for ConnectivityArrays
 *  that have variable ID types. First appends the initial_values then sets
 *  for half of n_IDs sets the values one at a time then sets the rest all at
 *  once.
 */
template< IndexType STRIDE, CellTypes TYPE >
void set( ConnectivityArray< STRIDE, TYPE >& connec, IndexType n_IDs,
          const IndexType* initial_values, const IndexType* values, 
          const IndexType* offsets=AXOM_NULLPTR, 
          const CellTypes* types=AXOM_NULLPTR )
{
  const IndexType initial_n_IDs = connec.getNumberOfIDs();
  
  /* Append the initial values. */
  append( connec, n_IDs, initial_values, offsets, types );

  IndexType half_n_IDs = n_IDs / 2;
  for ( IndexType ID = 0; ID < half_n_IDs; ++ID )
  {
    const IndexType cur_n_values = get_n_values< STRIDE >( ID, offsets );
    const IndexType* cur_values = get_value_ptr< STRIDE >( ID, values, offsets );
    const CellTypes type = get_type< TYPE >( ID, types );

    connec.set( cur_values, initial_n_IDs + ID );
    EXPECT_EQ( connec.getIDType( initial_n_IDs + ID ), type );
    EXPECT_EQ( connec.getNumberOfValuesForID( initial_n_IDs + ID ), cur_n_values );
    check_pointers( cur_values, connec[ initial_n_IDs + ID ], cur_n_values, 
                    __LINE__ );
  }

  const IndexType cur_ID = half_n_IDs;
  const IndexType* cur_values = get_value_ptr< STRIDE >( cur_ID, values, offsets );
  connec.setM( cur_values, initial_n_IDs + cur_ID, n_IDs - cur_ID );
  
  for ( IndexType ID = 0; ID < n_IDs; ++ID )
  {
    const IndexType cur_n_values = get_n_values< STRIDE >( ID, offsets);
    const IndexType* cur_values = get_value_ptr< STRIDE >( ID, values, offsets );
    const CellTypes cur_type = get_type< TYPE >( ID, types );
    EXPECT_EQ( connec.getIDType( initial_n_IDs + ID ), cur_type );
    EXPECT_EQ( connec.getNumberOfValuesForID( initial_n_IDs + ID ), cur_n_values );
    check_pointers( cur_values, connec[ initial_n_IDs + ID ], cur_n_values, 
                    __LINE__ );
  }
}

/*!
 * \brief Check that inserting the values into the given ConnectivityArray 
 *  functions as expected.
 *
 * \param [in/out] connec the ConnectivityArray to insert into.
 * \param [in] n_IDs the number of IDs to insert.
 * \param [in] values the array of values to set.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 *
 * \pre connec.empty()
 * \pre n_IDs % 6 == 0
 *
 * \note offsets only need to be specified for ConnectivityArrays that need 
 *  indirection and types only needs to be specified for ConnectivityArrays
 *  that have variable ID types. For n_IDs / 6 repetitions inserts an ID
 *  to the front, inside, and rear of the ConnectivityArray. Finally inserts
 *  n_IDs / 6 IDs each to the front, interior, and rear of the
 *  ConnectivityArray in one go.
 */
template< IndexType STRIDE, CellTypes TYPE >
void insert( ConnectivityArray< STRIDE, TYPE >& connec, IndexType n_IDs,
             const IndexType* values, const IndexType* offsets=AXOM_NULLPTR, 
             const CellTypes* types=AXOM_NULLPTR )
{
  SLIC_ERROR_IF( !connec.empty(), 
                 "Insertion test requires an empty ConnectivityArray." );
  SLIC_ERROR_IF( n_IDs % 6 != 0, 
                 "Insertion test requires n_IDs to be a multiple of 6." );

  const IndexType half_n_IDs = n_IDs / 2;
  const IndexType third_n_IDs = n_IDs / 3;
  const IndexType sixth_n_IDs = n_IDs / 6;

  /* Array of IDs the next ID to insert at the front, middle, and back. */
  IndexType IDs[3] = { third_n_IDs - 1, third_n_IDs, 2 * third_n_IDs };

  /* Array of positions at which to insert for the front, middle, and back. */
  IndexType insert_positions[3] = { 0, 1, 2 };
  
  for ( IndexType round = 0; round < sixth_n_IDs; ++round )
  {
    for ( IndexType i = 0; i < 3; ++i )
    {
      const IndexType cur_ID = IDs[ i ];
      const IndexType insert_pos = insert_positions[ i ];
      const IndexType cur_n_values = get_n_values< STRIDE >( cur_ID, offsets );
      const IndexType* cur_values = get_value_ptr< STRIDE >( cur_ID, values, offsets );
      const CellTypes type = get_type< TYPE >( cur_ID, types );

      connec.insert( cur_values, insert_pos, cur_n_values, type );
      EXPECT_EQ( connec.getNumberOfIDs(), 3 * round + i + 1 );
      EXPECT_EQ( connec.getIDType( insert_pos ), type );
      EXPECT_EQ( connec.getNumberOfValuesForID( insert_pos ), cur_n_values );
      check_pointers( cur_values, connec[ insert_pos ], cur_n_values, __LINE__ );
    }

    /* The middle insertion increases by two, the back by three. */
    insert_positions[1] += 2;
    insert_positions[2] += 3;

    /* The front ID decreases by one while the middle and rear increase by one. */
    IDs[0]--;
    IDs[1]++;
    IDs[2]++;
  }

  /* Insert the rest of the front values. */
  const IndexType front_ID = 0;
  const IndexType* front_offsets = get_offset_ptr( front_ID, offsets );
  const IndexType* front_values = get_value_ptr< STRIDE >( front_ID, values, offsets );
  const CellTypes* front_types = get_type_ptr( front_ID, types );
  connec.insertM( front_values , front_ID, sixth_n_IDs, front_offsets, 
                  front_types );
  EXPECT_EQ( connec.getNumberOfIDs(), 2 * third_n_IDs );

  /* Insert the rest of the middle values. */
  const IndexType middle_ID = half_n_IDs;
  const IndexType* middle_offsets = get_offset_ptr( middle_ID, offsets );
  const IndexType* middle_values = get_value_ptr< STRIDE >( middle_ID, values, offsets );
  const CellTypes* middle_types = get_type_ptr( middle_ID, types );
  connec.insertM( middle_values, middle_ID, sixth_n_IDs, middle_offsets, 
                  middle_types );
  EXPECT_EQ( connec.getNumberOfIDs(), 5 * sixth_n_IDs );

  /* Insert the rest of the back values. */
  const IndexType back_ID = 5 * sixth_n_IDs;
  const IndexType* back_offsets = get_offset_ptr( back_ID, offsets );
  const IndexType* back_values = get_value_ptr< STRIDE >( back_ID, values, offsets );
  const CellTypes* back_types = get_type_ptr( back_ID, types );
  connec.insertM( back_values, back_ID, sixth_n_IDs, back_offsets, back_types );
  EXPECT_EQ( connec.getNumberOfIDs(), n_IDs );

  /* Check that the values were correctly inserted. */
  for ( IndexType ID = 0; ID < n_IDs; ++ID )
  {
    const IndexType cur_n_values = get_n_values< STRIDE >( ID, offsets);
    const IndexType* cur_values = get_value_ptr< STRIDE >( ID, values, offsets );
    const CellTypes cur_type = get_type< TYPE >( ID, types );
    EXPECT_EQ( connec.getIDType( ID ), cur_type );
    EXPECT_EQ( connec.getNumberOfValuesForID( ID ), cur_n_values );
    check_pointers( cur_values, connec[ ID ], cur_n_values, __LINE__ );
  }
}

/*!
 * \brief Check that capacity of the given ConnectivityArray functions as 
 *  expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to.
 * \param [in] n_IDs the number of IDs to append.
 * \param [in] values the array of values to append.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 *
 * \pre connec.empty()
 *
 * \note offsets only need to be specified for ConnectivityArrays that need 
 *  indirection and types only needs to be specified for ConnectivityArrays
 *  that have variable ID types. First reserves space for half of n_IDs and half
 *  the values, then appends that many IDs one at a time. Then appends one extra
 *  ID and checks that the resizing occured as expected. Then shrinks the
 *  ConnectivityArray and appends the remaining IDs all at once.
 */
template< IndexType STRIDE, CellTypes TYPE >
void append_capacity( ConnectivityArray< STRIDE, TYPE >& connec, 
                      IndexType n_IDs, const IndexType* values, 
                      const IndexType* offsets=AXOM_NULLPTR, 
                      const CellTypes* types=AXOM_NULLPTR )
{
  SLIC_ERROR_IF( !connec.empty(), 
                 "Insertion test requires an empty ConnectivityArray." );

  const IndexType half_n_IDs = n_IDs / 2;
  const IndexType n_values = get_n_values< STRIDE >( 0, n_IDs, offsets );
  const IndexType first_half_n_values = get_n_values< STRIDE >( 0, half_n_IDs, offsets );
  IndexType cur_ID_size = connec.getNumberOfIDs();
  IndexType cur_value_size = connec.getNumberOfValues();

  EXPECT_EQ( cur_ID_size, 0 );
  EXPECT_EQ( cur_value_size, 0 );

  connec.reserve( half_n_IDs, first_half_n_values );
  const IndexType* cur_values_ptr = connec.getValuePtr();
  const IndexType* cur_offsets_ptr = connec.getOffsetPtr();
  const CellTypes* cur_types_ptr = connec.getTypePtr();
  IndexType cur_ID_capacity = connec.getIDCapacity();
  IndexType cur_value_capacity = connec.getValueCapacity();

  EXPECT_NE( cur_values_ptr, AXOM_NULLPTR );
  EXPECT_EQ( cur_ID_capacity, half_n_IDs );
  EXPECT_EQ( cur_value_capacity, first_half_n_values );

  /* Append the first half of the IDs, no resize should occur. */
  for ( IndexType ID = 0; ID < half_n_IDs; ++ID )
  {
    const IndexType cur_n_values = get_n_values< STRIDE >( ID, offsets );
    const IndexType* cur_values = get_value_ptr< STRIDE >( ID, values, offsets );
    const CellTypes type = get_type< TYPE >( ID, types );

    connec.append( cur_values, cur_n_values, type );
    cur_ID_size += 1;
    cur_value_size += cur_n_values;
    EXPECT_EQ( connec.getIDCapacity(), cur_ID_capacity );
    EXPECT_EQ( connec.getValueCapacity(), cur_value_capacity );
    EXPECT_EQ( connec.getValuePtr(), cur_values_ptr );
    EXPECT_EQ( connec.getOffsetPtr(), cur_offsets_ptr );
    EXPECT_EQ( connec.getTypePtr(), cur_types_ptr );
    EXPECT_EQ( connec.getNumberOfIDs(), cur_ID_size );
    EXPECT_EQ( connec.getNumberOfValues(), cur_value_size );
    EXPECT_EQ( connec.getIDType( ID ), type );
    EXPECT_EQ( connec.getNumberOfValuesForID( ID ), cur_n_values );
    check_pointers( cur_values, connec[ ID ], cur_n_values, __LINE__ );
  }

  /* Append one more value, should trigger a resize. */
  IndexType cur_ID = half_n_IDs;
  IndexType cur_n_values = get_n_values< STRIDE >( cur_ID, offsets );
  const IndexType* cur_values = get_value_ptr< STRIDE >( cur_ID, values, offsets );
  CellTypes type = get_type< TYPE >( cur_ID, types );
  
  cur_ID_capacity = calc_ID_capacity( connec, 1 );
  cur_value_capacity = calc_value_capacity( connec, cur_n_values );
  connec.append( cur_values, cur_n_values, type );
  cur_ID_size += 1;
  cur_value_size += cur_n_values;

  EXPECT_EQ( connec.getIDCapacity(), cur_ID_capacity );
  EXPECT_EQ( connec.getValueCapacity(), cur_value_capacity );
  EXPECT_EQ( connec.getNumberOfIDs(), cur_ID_size );
  EXPECT_EQ( connec.getNumberOfValues(), cur_value_size );
  EXPECT_EQ( connec.getIDType( cur_ID ), type );
  EXPECT_EQ( connec.getNumberOfValuesForID( cur_ID ), cur_n_values );
  check_pointers( cur_values, connec[ cur_ID ], cur_n_values, __LINE__ );

  /* Shrink. */
  connec.shrink();
  cur_ID_capacity = cur_ID_size;
  cur_value_capacity = cur_value_size;
  EXPECT_EQ( connec.getIDCapacity(), cur_ID_capacity );
  EXPECT_EQ( connec.getValueCapacity(), cur_value_capacity );

  /* Append the rest of the values all at once, should trigger a resize. */
  cur_ID = half_n_IDs + 1;
  cur_values = get_value_ptr< STRIDE >( cur_ID, values, offsets );
  cur_n_values = get_n_values< STRIDE >( cur_ID, n_IDs, offsets );
  const IndexType* cur_offsets = get_offset_ptr( cur_ID, offsets );
  const CellTypes* cur_types = get_type_ptr( cur_ID, types );
  
  
  cur_ID_capacity = calc_ID_capacity( connec, n_IDs - cur_ID );
  cur_value_capacity = calc_value_capacity( connec, cur_n_values );
  connec.appendM( cur_values, n_IDs - cur_ID, cur_offsets, cur_types );
  cur_ID_size += n_IDs - cur_ID;
  cur_value_size += cur_n_values;

  EXPECT_EQ( connec.getIDCapacity(), cur_ID_capacity );
  EXPECT_EQ( connec.getValueCapacity(), cur_value_capacity );
  EXPECT_EQ( connec.getNumberOfIDs(), cur_ID_size );
  EXPECT_EQ( connec.getNumberOfValues(), cur_value_size );
  EXPECT_EQ( connec.getNumberOfIDs(), n_IDs );
  EXPECT_EQ( connec.getNumberOfValues(), n_values );

  for ( IndexType ID = 0; ID < n_IDs; ++ID )
  {
    cur_n_values = get_n_values< STRIDE >( ID, offsets);
    cur_values = get_value_ptr< STRIDE >( ID, values, offsets );
    type = get_type< TYPE >( ID, types );
    EXPECT_EQ( connec.getIDType( ID ), type );
    check_pointers( cur_values, connec[ ID ], cur_n_values, __LINE__ );
  }

  check_pointers( connec.getValuePtr(), values, n_values, __LINE__ );
  check_pointers( connec.getOffsetPtr(), offsets, n_IDs + 1, __LINE__ );
  check_pointers( connec.getTypePtr(), types, n_IDs, __LINE__ );
}


}   /* end namespace internal */


/*******************************************************************************
 *                          Append tests                                       *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, NoIndirectionNativeAppend )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }
  
  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray< vertex_stride, vertex > native_vertex;
  EXPECT_FALSE( native_vertex.isExternal() );
  EXPECT_FALSE( native_vertex.isInSidre() );
  EXPECT_TRUE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< hex_stride, CellTypes::HEX > native_hex;
  EXPECT_FALSE( native_hex.isExternal() );
  EXPECT_FALSE( native_hex.isInSidre() );
  EXPECT_TRUE( native_hex.empty() );
  EXPECT_EQ( native_hex.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( native_hex.getGroup(), AXOM_NULLPTR );
#endif

  /* Append the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::append( native_vertex, n_IDs, values );
    internal::append( native_hex, n_IDs, values );
  }

  /* Check that the values were appended properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::check_pointers( values, native_vertex[ cur_offset ], 
                              n_IDs * vertex_stride, __LINE__ );
    internal::check_pointers( values, native_hex[ cur_offset ], 
                              n_IDs * hex_stride, __LINE__ );
    cur_offset += n_IDs;
  }

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_vertex.getNumberOfValues(), total_IDs * vertex_stride);

  EXPECT_FALSE( native_hex.empty() );
  EXPECT_EQ( native_hex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_hex.getNumberOfValues(), total_IDs * hex_stride);

  delete[] values;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, NoIndirectionExternalAppend )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[ total_IDs * vertex_stride ];
  ConnectivityArray< vertex_stride, vertex > ext_vertex( 0, ext_vertex_values, 
                                                         total_IDs );
  EXPECT_TRUE( ext_vertex.isExternal() );
  EXPECT_FALSE( ext_vertex.isInSidre() );
  EXPECT_TRUE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), 0);
  EXPECT_EQ( ext_vertex.getValuePtr(), ext_vertex_values );

  IndexType* ext_hex_values = new IndexType[ total_IDs * hex_stride ];
  ConnectivityArray< hex_stride, hex > ext_hex( 0, ext_hex_values, 
                                                total_IDs );
  EXPECT_TRUE( ext_hex.isExternal() );
  EXPECT_FALSE( ext_hex.isInSidre() );
  EXPECT_TRUE( ext_hex.empty() );
  EXPECT_EQ( ext_hex.getNumberOfIDs(), 0);
  EXPECT_EQ( ext_hex.getValuePtr(), ext_hex_values );

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( ext_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( ext_hex.getGroup(), AXOM_NULLPTR );
#endif

  /* Append the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::append( ext_vertex, n_IDs, values );
    internal::append( ext_hex, n_IDs, values );
  }

  /* Check that the values were appended properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::check_pointers( values, ext_vertex[ cur_offset ], 
                              n_IDs * vertex_stride, __LINE__ );
    internal::check_pointers( values, ext_hex[ cur_offset ], 
                              n_IDs * hex_stride, __LINE__ );
    cur_offset += n_IDs;
  }

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_vertex.getNumberOfValues(), total_IDs * vertex_stride);
  EXPECT_EQ( ext_vertex.getValuePtr(), ext_vertex_values );

  EXPECT_FALSE( ext_hex.empty() );
  EXPECT_EQ( ext_hex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_hex.getNumberOfValues(), total_IDs * hex_stride);
  EXPECT_EQ( ext_hex.getValuePtr(), ext_hex_values );


  /* Check that the external ConnectivityArrays cannot append any more. */
  EXPECT_DEATH_IF_SUPPORTED( ext_vertex.append( values ), IGNORE_OUTPUT ); 
  EXPECT_DEATH_IF_SUPPORTED( ext_hex.append( values ), IGNORE_OUTPUT ); 

  /* Check that the external constructor functions properly. */
  ConnectivityArray< vertex_stride, vertex > 
              ext_vertex_cpy( ext_vertex.getNumberOfIDs(), ext_vertex_values );
  internal::check_equality( ext_vertex, ext_vertex_cpy, __LINE__ );

  ConnectivityArray< hex_stride, hex > 
                        ext_hex_cpy( ext_hex.getNumberOfIDs(), ext_hex_values );
  internal::check_equality( ext_hex, ext_hex_cpy, __LINE__ );

  delete[] values;
  delete[] ext_vertex_values;
  delete[] ext_hex_values;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, NoIndirectionSidreAppend )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }
  
  /* Create the sidre storage ConnectivityArrays to be tested. */
  ConnectivityArray< vertex_stride, vertex > 
                          sidre_vertex( root->createGroup( "vertex" ), "test" );
  EXPECT_TRUE( sidre_vertex.isInSidre() );
  EXPECT_FALSE( sidre_vertex.isExternal() );
  EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
  EXPECT_TRUE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< hex_stride, hex >
                          sidre_hex( root->createGroup( "hex" ), "test" );
  EXPECT_TRUE( sidre_hex.isInSidre() );
  EXPECT_FALSE( sidre_hex.isExternal() );
  EXPECT_EQ( sidre_hex.getGroup(), root->getGroup( "hex") );
  EXPECT_TRUE( sidre_hex.empty() );
  EXPECT_EQ( sidre_hex.getNumberOfIDs(), 0);

  /* Append the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::append( sidre_vertex, n_IDs, values );
    internal::append( sidre_hex, n_IDs, values );
  }

  /* Check that the values were appended properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::check_pointers( values, sidre_vertex[ cur_offset ], 
                              n_IDs * vertex_stride, __LINE__ );
    internal::check_pointers( values, sidre_hex[ cur_offset ], 
                              n_IDs * hex_stride, __LINE__ );
    cur_offset += n_IDs;
  }

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_vertex.getNumberOfValues(), total_IDs * vertex_stride);

  EXPECT_FALSE( sidre_hex.empty() );
  EXPECT_EQ( sidre_hex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_hex.getNumberOfValues(), total_IDs * hex_stride);

  /* Check that restoring from sidre functions properly */
  sidre::Group* vertex_group = 
                        const_cast< sidre::Group* >( sidre_vertex.getGroup() );
  ConnectivityArray< vertex_stride, vertex > sidre_vertex_cpy( vertex_group );
  internal::check_equality( sidre_vertex, sidre_vertex_cpy, __LINE__ );

  sidre::Group* hex_group = 
                        const_cast< sidre::Group* >( sidre_hex.getGroup() );
  ConnectivityArray< hex_stride, hex > sidre_hex_cpy( hex_group );
  internal::check_equality( sidre_hex, sidre_hex_cpy, __LINE__ );

  delete[] values;

#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, IndirectionNativeAppend )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * max_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ max_IDs + 1 ];
  CellTypes* types = new CellTypes[ max_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < max_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  /* Calculate the total number of IDs and values. */
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;
  IndexType total_values = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    total_values += (n_IDs / 2) * (vertex_stride + hex_stride);
    if (n_IDs % 2 == 1)
    {
      total_values += vertex_stride;
    }
  }
  
  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray< NEED_INDIRECTION, vertex > native_vertex;
  EXPECT_FALSE( native_vertex.isExternal() );
  EXPECT_FALSE( native_vertex.isInSidre() );
  EXPECT_TRUE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< NEED_INDIRECTION, mixed > native_mixed;
  EXPECT_FALSE( native_mixed.isExternal() );
  EXPECT_FALSE( native_mixed.isInSidre() );
  EXPECT_TRUE( native_mixed.empty() );
  EXPECT_EQ( native_mixed.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( native_mixed.getGroup(), AXOM_NULLPTR );
#endif

  /* Append the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::append( native_vertex, n_IDs, values, offsets );
    internal::append( native_mixed, n_IDs, values, offsets, types );
  }

  /* Check that the values were appended properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    for ( IndexType ID = 0; ID < n_IDs; ++ID )
    {
      const IndexType* cur_values = values + offsets[ ID ];
      const IndexType n_values = offsets[ ID + 1 ] - offsets[ ID ];
      internal::check_pointers( cur_values, native_vertex[ cur_offset + ID ], 
                                n_values, __LINE__ );
      internal::check_pointers( cur_values, native_mixed[ cur_offset + ID ], 
                                n_values, __LINE__ );
    }

    cur_offset += n_IDs;
  }

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_vertex.getNumberOfValues(), total_values);

  EXPECT_FALSE( native_mixed.empty() );
  EXPECT_EQ( native_mixed.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_mixed.getNumberOfValues(), total_values);

  delete[] values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, IndirectionExternalAppend )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * max_IDs / 2;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ max_IDs + 1 ];
  CellTypes* types = new CellTypes[ max_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < max_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  /* Calculate the total number of IDs and values. */
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;
  IndexType total_values = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    total_values += (n_IDs / 2) * (vertex_stride + hex_stride);
    if (n_IDs % 2 == 1)
    {
      total_values += vertex_stride;
    }
  }
  
  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[ total_values ];
  IndexType* ext_vertex_offsets = new IndexType[ total_IDs + 1 ];
  ConnectivityArray< NEED_INDIRECTION, vertex > 
                          ext_vertex( 0, ext_vertex_values, ext_vertex_offsets,
                                      total_IDs, total_values );
  EXPECT_TRUE( ext_vertex.isExternal() );
  EXPECT_FALSE( ext_vertex.isInSidre() );
  EXPECT_TRUE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), 0);
  EXPECT_EQ( ext_vertex.getValuePtr(), ext_vertex_values );
  EXPECT_EQ( ext_vertex.getOffsetPtr(), ext_vertex_offsets );

  IndexType* ext_mixed_values = new IndexType[ total_IDs * hex_stride ];
  IndexType* ext_mixed_offsets = new IndexType[ total_IDs + 1 ];
  CellTypes* ext_mixed_types = new CellTypes[ total_IDs ];
  ConnectivityArray< NEED_INDIRECTION, mixed > 
                          ext_mixed( 0, ext_mixed_values, ext_mixed_offsets, 
                                     ext_mixed_types, total_IDs, total_values );
  EXPECT_TRUE( ext_mixed.isExternal() );
  EXPECT_FALSE( ext_mixed.isInSidre() );
  EXPECT_TRUE( ext_mixed.empty() );
  EXPECT_EQ( ext_mixed.getNumberOfIDs(), 0);
  EXPECT_EQ( ext_mixed.getValuePtr(), ext_mixed_values );
  EXPECT_EQ( ext_mixed.getOffsetPtr(), ext_mixed_offsets );
  EXPECT_EQ( ext_mixed.getTypePtr(), ext_mixed_types );

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( ext_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( ext_mixed.getGroup(), AXOM_NULLPTR );
#endif

  /* Append the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::append( ext_vertex, n_IDs, values, offsets );
    internal::append( ext_mixed, n_IDs, values, offsets, types );
  }

  /* Check that the values were appended properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    for ( IndexType ID = 0; ID < n_IDs; ++ID )
    {
      const IndexType* cur_values = values + offsets[ ID ];
      const IndexType n_values = offsets[ ID + 1 ] - offsets[ ID ];
      internal::check_pointers( cur_values, ext_vertex[ cur_offset + ID ], 
                                n_values, __LINE__ );
      internal::check_pointers( cur_values, ext_mixed[ cur_offset + ID ], 
                                n_values, __LINE__ );
    }

    cur_offset += n_IDs;
  }

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_vertex.getNumberOfValues(), total_values);
  EXPECT_EQ( ext_vertex.getValuePtr(), ext_vertex_values );
  EXPECT_EQ( ext_vertex.getOffsetPtr(), ext_vertex_offsets );

  EXPECT_FALSE( ext_mixed.empty() );
  EXPECT_EQ( ext_mixed.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_mixed.getNumberOfValues(), total_values);
  EXPECT_EQ( ext_mixed.getValuePtr(), ext_mixed_values );
  EXPECT_EQ( ext_mixed.getOffsetPtr(), ext_mixed_offsets );
  EXPECT_EQ( ext_mixed.getTypePtr(), ext_mixed_types );

  /* Check that the external ConnectivityArrays cannot append any more. */
  EXPECT_DEATH_IF_SUPPORTED( ext_vertex.append( values, 1 ), IGNORE_OUTPUT ); 
  EXPECT_DEATH_IF_SUPPORTED( ext_mixed.append( values,  1, vertex ), IGNORE_OUTPUT );

  /* Check that the external constructor functions properly. */
  ConnectivityArray< NEED_INDIRECTION, vertex > 
              ext_vertex_cpy( ext_vertex.getNumberOfIDs(), ext_vertex_values, 
                              ext_vertex_offsets );
  internal::check_equality( ext_vertex, ext_vertex_cpy, __LINE__ );

  ConnectivityArray< NEED_INDIRECTION, mixed > 
                        ext_mixed_cpy( ext_mixed.getNumberOfIDs(), ext_mixed_values, 
                                       ext_mixed_offsets, ext_mixed_types );
  internal::check_equality( ext_mixed, ext_mixed_cpy, __LINE__ );

  delete[] values;
  delete[] offsets;
  delete[] types;
  delete[] ext_vertex_values;
  delete[] ext_vertex_offsets;
  delete[] ext_mixed_values;
  delete[] ext_mixed_offsets;
  delete[] ext_mixed_types;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, IndirectionSidreAppend )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * max_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ max_IDs + 1 ];
  CellTypes* types = new CellTypes[ max_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < max_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  /* Calculate the total number of IDs and values. */
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;
  IndexType total_values = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    total_values += (n_IDs / 2) * (vertex_stride + hex_stride);
    if (n_IDs % 2 == 1)
    {
      total_values += vertex_stride;
    }
  }

  /* Create the sidre storage ConnectivityArrays to be tested. */
  ConnectivityArray< NEED_INDIRECTION, vertex > 
                          sidre_vertex( root->createGroup( "vertex" ), "test" );
  EXPECT_TRUE( sidre_vertex.isInSidre() );
  EXPECT_FALSE( sidre_vertex.isExternal() );
  EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
  EXPECT_TRUE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< NEED_INDIRECTION, mixed >
                          sidre_mixed( root->createGroup( "mixed" ), "test" );
  EXPECT_TRUE( sidre_mixed.isInSidre() );
  EXPECT_FALSE( sidre_mixed.isExternal() );
  EXPECT_EQ( sidre_mixed.getGroup(), root->getGroup( "mixed" ) );
  EXPECT_TRUE( sidre_mixed.empty() );
  EXPECT_EQ( sidre_mixed.getNumberOfIDs(), 0);

  /* Append the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::append( sidre_vertex, n_IDs, values, offsets );
    internal::append( sidre_mixed, n_IDs, values, offsets, types );
  }

  /* Check that the values were appended properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    for ( IndexType ID = 0; ID < n_IDs; ++ID )
    {
      const IndexType* cur_values = values + offsets[ ID ];
      const IndexType n_values = offsets[ ID + 1 ] - offsets[ ID ];
      internal::check_pointers( cur_values, sidre_vertex[ cur_offset + ID ], 
                                n_values, __LINE__ );
      internal::check_pointers( cur_values, sidre_mixed[ cur_offset + ID ], 
                                n_values, __LINE__ );
    }

    cur_offset += n_IDs;
  }

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_vertex.getNumberOfValues(), total_values);

  EXPECT_FALSE( sidre_mixed.empty() );
  EXPECT_EQ( sidre_mixed.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_mixed.getNumberOfValues(), total_values);

  /* Check that restoring from sidre functions properly */
  sidre::Group* vertex_group = 
                        const_cast< sidre::Group* >( sidre_vertex.getGroup() );
  ConnectivityArray< NEED_INDIRECTION, vertex > sidre_vertex_cpy( vertex_group );
  internal::check_equality( sidre_vertex, sidre_vertex_cpy, __LINE__ );

  sidre::Group* mixed_group = 
                        const_cast< sidre::Group* >( sidre_mixed.getGroup() );
  ConnectivityArray< NEED_INDIRECTION, mixed > sidre_mixed_cpy( mixed_group );
  internal::check_equality( sidre_mixed, sidre_mixed_cpy, __LINE__ );

  delete[] values;
  delete[] offsets;
  delete[] types;

#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}


/*******************************************************************************
 *                          Set tests                                          *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, NoIndirectionNativeSet )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;

  /* Allocate and populate the values buffer. */
  IndexType* initial_values = new IndexType[ max_values ];
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    initial_values[ i ] = -i;
    values[ i ] = i;
  }
  
  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray< vertex_stride, vertex > native_vertex;
  EXPECT_FALSE( native_vertex.isExternal() );
  EXPECT_FALSE( native_vertex.isInSidre() );
  EXPECT_TRUE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< hex_stride, CellTypes::HEX > native_hex;
  EXPECT_FALSE( native_hex.isExternal() );
  EXPECT_FALSE( native_hex.isInSidre() );
  EXPECT_TRUE( native_hex.empty() );
  EXPECT_EQ( native_hex.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( native_hex.getGroup(), AXOM_NULLPTR );
#endif

  /* Append and set the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::set( native_vertex, n_IDs, initial_values, values );
    internal::set( native_hex, n_IDs, initial_values, values );
  }

  /* Check that the values were appended and set properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::check_pointers( values, native_vertex[ cur_offset ], 
                              n_IDs * vertex_stride, __LINE__ );
    internal::check_pointers( values, native_hex[ cur_offset ], 
                              n_IDs * hex_stride, __LINE__ );
    cur_offset += n_IDs;
  }

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_vertex.getNumberOfValues(), total_IDs * vertex_stride);

  EXPECT_FALSE( native_hex.empty() );
  EXPECT_EQ( native_hex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_hex.getNumberOfValues(), total_IDs * hex_stride);

  delete[] values;
  delete[] initial_values;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, NoIndirectionExternalSet )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  IndexType* initial_values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
    initial_values[ i ] = -i;
  }

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[ total_IDs * vertex_stride ];
  ConnectivityArray< vertex_stride, vertex > ext_vertex( 0, ext_vertex_values, 
                                                         total_IDs );
  EXPECT_TRUE( ext_vertex.isExternal() );
  EXPECT_FALSE( ext_vertex.isInSidre() );
  EXPECT_TRUE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), 0);

  IndexType* ext_hex_values = new IndexType[ total_IDs * hex_stride ];
  ConnectivityArray< hex_stride, hex > ext_hex( 0, ext_hex_values, 
                                                total_IDs );
  EXPECT_TRUE( ext_hex.isExternal() );
  EXPECT_FALSE( ext_hex.isInSidre() );
  EXPECT_TRUE( ext_hex.empty() );
  EXPECT_EQ( ext_hex.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( ext_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( ext_hex.getGroup(), AXOM_NULLPTR );
#endif

  /* Append and set the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::set( ext_vertex, n_IDs, initial_values, values );
    internal::set( ext_hex, n_IDs, initial_values, values );
  }

  /* Check that the values were appended and set properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::check_pointers( values, ext_vertex[ cur_offset ], 
                              n_IDs * vertex_stride, __LINE__ );
    internal::check_pointers( values, ext_hex[ cur_offset ], 
                              n_IDs * hex_stride, __LINE__ );
    cur_offset += n_IDs;
  }

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_vertex.getNumberOfValues(), total_IDs * vertex_stride);

  EXPECT_FALSE( ext_hex.empty() );
  EXPECT_EQ( ext_hex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_hex.getNumberOfValues(), total_IDs * hex_stride);

  delete[] values;
  delete[] initial_values;
  delete[] ext_vertex_values;
  delete[] ext_hex_values;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, NoIndirectionSidreSet )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  IndexType* initial_values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
    initial_values[ i ] = -i;
  }
  
  /* Create the sidre storage ConnectivityArrays to be tested. */
  ConnectivityArray< vertex_stride, vertex > 
                          sidre_vertex( root->createGroup( "vertex" ), "test" );
  EXPECT_TRUE( sidre_vertex.isInSidre() );
  EXPECT_FALSE( sidre_vertex.isExternal() );
  EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
  EXPECT_TRUE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< hex_stride, hex >
                          sidre_hex( root->createGroup( "hex" ), "test" );
  EXPECT_TRUE( sidre_hex.isInSidre() );
  EXPECT_FALSE( sidre_hex.isExternal() );
  EXPECT_EQ( sidre_hex.getGroup(), root->getGroup( "hex") );
  EXPECT_TRUE( sidre_hex.empty() );
  EXPECT_EQ( sidre_hex.getNumberOfIDs(), 0);

  /* Append and set the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::set( sidre_vertex, n_IDs, initial_values, values );
    internal::set( sidre_hex, n_IDs, initial_values, values );
  }

  /* Check that the values were appended and set properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::check_pointers( values, sidre_vertex[ cur_offset ], 
                              n_IDs * vertex_stride, __LINE__ );
    internal::check_pointers( values, sidre_hex[ cur_offset ], 
                              n_IDs * hex_stride, __LINE__ );
    cur_offset += n_IDs;
  }

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_vertex.getNumberOfValues(), total_IDs * vertex_stride);

  EXPECT_FALSE( sidre_hex.empty() );
  EXPECT_EQ( sidre_hex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_hex.getNumberOfValues(), total_IDs * hex_stride);

  /* Check that restoring from sidre functions properly */
  sidre::Group* vertex_group = 
                        const_cast< sidre::Group* >( sidre_vertex.getGroup() );
  ConnectivityArray< vertex_stride, vertex > sidre_vertex_cpy( vertex_group );
  internal::check_equality( sidre_vertex, sidre_vertex_cpy, __LINE__ );

  sidre::Group* hex_group = 
                        const_cast< sidre::Group* >( sidre_hex.getGroup() );
  ConnectivityArray< hex_stride, hex > sidre_hex_cpy( hex_group );
  internal::check_equality( sidre_hex, sidre_hex_cpy, __LINE__ );

  delete[] values;
  delete[] initial_values;

#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, IndirectionNativeSet )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * max_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  IndexType* initial_values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
    initial_values[ i ] = -i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ max_IDs + 1 ];
  CellTypes* types = new CellTypes[ max_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < max_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  /* Calculate the total number of IDs and values. */
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;
  IndexType total_values = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    total_values += (n_IDs / 2) * (vertex_stride + hex_stride);
    if (n_IDs % 2 == 1)
    {
      total_values += vertex_stride;
    }
  }
  
  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray< NEED_INDIRECTION, vertex > native_vertex;
  EXPECT_FALSE( native_vertex.isExternal() );
  EXPECT_FALSE( native_vertex.isInSidre() );
  EXPECT_TRUE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< NEED_INDIRECTION, mixed > native_mixed;
  EXPECT_FALSE( native_mixed.isExternal() );
  EXPECT_FALSE( native_mixed.isInSidre() );
  EXPECT_TRUE( native_mixed.empty() );
  EXPECT_EQ( native_mixed.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( native_mixed.getGroup(), AXOM_NULLPTR );
#endif

  /* Append and set the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::set( native_vertex, n_IDs, initial_values, values, offsets );
    internal::set( native_mixed, n_IDs, initial_values, values, offsets, types );
  }

  /* Check that the values were appended and set properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    for ( IndexType ID = 0; ID < n_IDs; ++ID )
    {
      const IndexType* cur_values = values + offsets[ ID ];
      const IndexType n_values = offsets[ ID + 1 ] - offsets[ ID ];
      internal::check_pointers( cur_values, native_vertex[ cur_offset + ID ], 
                                n_values, __LINE__ );
      internal::check_pointers( cur_values, native_mixed[ cur_offset + ID ], 
                                n_values, __LINE__ );
    }

    cur_offset += n_IDs;
  }

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_vertex.getNumberOfValues(), total_values);

  EXPECT_FALSE( native_mixed.empty() );
  EXPECT_EQ( native_mixed.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( native_mixed.getNumberOfValues(), total_values);

  delete[] values;
  delete[] initial_values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, IndirectionExternalSet )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * max_IDs / 2;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  IndexType* initial_values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
    initial_values[ i ] = -i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ max_IDs + 1 ];
  CellTypes* types = new CellTypes[ max_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < max_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  /* Calculate the total number of IDs and values. */
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;
  IndexType total_values = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    total_values += (n_IDs / 2) * (vertex_stride + hex_stride);
    if (n_IDs % 2 == 1)
    {
      total_values += vertex_stride;
    }
  }
  
  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[ total_values ];
  IndexType* ext_vertex_offsets = new IndexType[ total_IDs + 1 ];
  ConnectivityArray< NEED_INDIRECTION, vertex > 
                          ext_vertex( 0, ext_vertex_values, ext_vertex_offsets,
                                      total_IDs, total_values );
  EXPECT_TRUE( ext_vertex.isExternal() );
  EXPECT_FALSE( ext_vertex.isInSidre() );
  EXPECT_TRUE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), 0);

  IndexType* ext_mixed_values = new IndexType[ total_IDs * hex_stride ];
  IndexType* ext_mixed_offsets = new IndexType[ total_IDs + 1 ];
  CellTypes* ext_mixed_types = new CellTypes[ total_IDs ];
  ConnectivityArray< NEED_INDIRECTION, mixed > 
                          ext_mixed( 0, ext_mixed_values, ext_mixed_offsets, 
                                     ext_mixed_types, total_IDs, total_values );
  EXPECT_TRUE( ext_mixed.isExternal() );
  EXPECT_FALSE( ext_mixed.isInSidre() );
  EXPECT_TRUE( ext_mixed.empty() );
  EXPECT_EQ( ext_mixed.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( ext_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( ext_mixed.getGroup(), AXOM_NULLPTR );
#endif

  /* Append and set the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::set( ext_vertex, n_IDs, initial_values, values, offsets );
    internal::set( ext_mixed, n_IDs, initial_values, values, offsets, types );
  }

  /* Check that the values were appended and set properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    for ( IndexType ID = 0; ID < n_IDs; ++ID )
    {
      const IndexType* cur_values = values + offsets[ ID ];
      const IndexType n_values = offsets[ ID + 1 ] - offsets[ ID ];
      internal::check_pointers( cur_values, ext_vertex[ cur_offset + ID ], 
                                n_values, __LINE__ );
      internal::check_pointers( cur_values, ext_mixed[ cur_offset + ID ], 
                                n_values, __LINE__ );
    }

    cur_offset += n_IDs;
  }

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_vertex.getNumberOfValues(), total_values);

  EXPECT_FALSE( ext_mixed.empty() );
  EXPECT_EQ( ext_mixed.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( ext_mixed.getNumberOfValues(), total_values);

  delete[] values;
  delete[] offsets;
  delete[] types;
  delete[] ext_vertex_values;
  delete[] ext_vertex_offsets;
  delete[] ext_mixed_values;
  delete[] ext_mixed_offsets;
  delete[] ext_mixed_types;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, IndirectionSidreSet )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * max_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  IndexType* initial_values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
    initial_values[ i ] = -i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ max_IDs + 1 ];
  CellTypes* types = new CellTypes[ max_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < max_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  /* Calculate the total number of IDs and values. */
  constexpr IndexType total_IDs = ( max_IDs * ( max_IDs + 1 ) ) / 2;
  IndexType total_values = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    total_values += (n_IDs / 2) * (vertex_stride + hex_stride);
    if (n_IDs % 2 == 1)
    {
      total_values += vertex_stride;
    }
  }

  /* Create the sidre storage ConnectivityArrays to be tested. */
  ConnectivityArray< NEED_INDIRECTION, vertex > 
                          sidre_vertex( root->createGroup( "vertex" ), "test" );
  EXPECT_TRUE( sidre_vertex.isInSidre() );
  EXPECT_FALSE( sidre_vertex.isExternal() );
  EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
  EXPECT_TRUE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< NEED_INDIRECTION, mixed >
                          sidre_mixed( root->createGroup( "mixed" ), "test" );
  EXPECT_TRUE( sidre_mixed.isInSidre() );
  EXPECT_FALSE( sidre_mixed.isExternal() );
  EXPECT_EQ( sidre_mixed.getGroup(), root->getGroup( "mixed" ) );
  EXPECT_TRUE( sidre_mixed.empty() );
  EXPECT_EQ( sidre_mixed.getNumberOfIDs(), 0);



  /* Append and set the values */
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    internal::set( sidre_vertex, n_IDs, initial_values, values, offsets );
    internal::set( sidre_mixed, n_IDs, initial_values, values, offsets, types );
  }

  /* Check that the values were appended and set properly */
  IndexType cur_offset = 0;
  for ( IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs )
  {
    for ( IndexType ID = 0; ID < n_IDs; ++ID )
    {
      const IndexType* cur_values = values + offsets[ ID ];
      const IndexType n_values = offsets[ ID + 1 ] - offsets[ ID ];
      internal::check_pointers( cur_values, sidre_vertex[ cur_offset + ID ], 
                                n_values, __LINE__ );
      internal::check_pointers( cur_values, sidre_mixed[ cur_offset + ID ], 
                                n_values, __LINE__ );
    }

    cur_offset += n_IDs;
  }

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_vertex.getNumberOfValues(), total_values);

  EXPECT_FALSE( sidre_mixed.empty() );
  EXPECT_EQ( sidre_mixed.getNumberOfIDs(), total_IDs);
  EXPECT_EQ( sidre_mixed.getNumberOfValues(), total_values);

  /* Check that restoring from sidre functions properly */
  sidre::Group* vertex_group = 
                        const_cast< sidre::Group* >( sidre_vertex.getGroup() );
  ConnectivityArray< NEED_INDIRECTION, vertex > sidre_vertex_cpy( vertex_group );
  internal::check_equality( sidre_vertex, sidre_vertex_cpy, __LINE__ );

  sidre::Group* mixed_group = 
                        const_cast< sidre::Group* >( sidre_mixed.getGroup() );
  ConnectivityArray< NEED_INDIRECTION, mixed > sidre_mixed_cpy( mixed_group );
  internal::check_equality( sidre_mixed, sidre_mixed_cpy, __LINE__ );

  delete[] values;
  delete[] initial_values;
  delete[] offsets;
  delete[] types;

#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}


/*******************************************************************************
 *                          Insert tests                                       *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, NoIndirectionNativeInsert )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType n_IDs = 500 * 6;
  constexpr IndexType n_values = hex_stride * n_IDs;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ n_values ];
  for ( IndexType i = 0; i < n_values; ++i )
  {
    values[ i ] = i;
  }
  
  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray< vertex_stride, vertex > native_vertex;
  EXPECT_FALSE( native_vertex.isExternal() );
  EXPECT_FALSE( native_vertex.isInSidre() );
  EXPECT_TRUE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< hex_stride, CellTypes::HEX > native_hex;
  EXPECT_FALSE( native_hex.isExternal() );
  EXPECT_FALSE( native_hex.isInSidre() );
  EXPECT_TRUE( native_hex.empty() );
  EXPECT_EQ( native_hex.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( native_hex.getGroup(), AXOM_NULLPTR );
#endif

  /* Insert the values */
  internal::insert( native_vertex, n_IDs, values );
  internal::insert( native_hex, n_IDs, values );

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( native_vertex.getNumberOfValues(), n_IDs * vertex_stride);

  EXPECT_FALSE( native_hex.empty() );
  EXPECT_EQ( native_hex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( native_hex.getNumberOfValues(), n_IDs * hex_stride);

  delete[] values;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, NoIndirectionExternalInsert )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType n_IDs = 500 * 6;
  constexpr IndexType n_values = hex_stride * n_IDs;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ n_values ];
  for ( IndexType i = 0; i < n_values; ++i )
  {
    values[ i ] = i;
  }

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[ n_IDs * vertex_stride ];
  ConnectivityArray< vertex_stride, vertex > ext_vertex( 0, ext_vertex_values, 
                                                         n_IDs );
  EXPECT_TRUE( ext_vertex.isExternal() );
  EXPECT_FALSE( ext_vertex.isInSidre() );
  EXPECT_TRUE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), 0);

  IndexType* ext_hex_values = new IndexType[ n_IDs * hex_stride ];
  ConnectivityArray< hex_stride, hex > ext_hex( 0, ext_hex_values, 
                                                n_IDs );
  EXPECT_TRUE( ext_hex.isExternal() );
  EXPECT_FALSE( ext_hex.isInSidre() );
  EXPECT_TRUE( ext_hex.empty() );
  EXPECT_EQ( ext_hex.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( ext_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( ext_hex.getGroup(), AXOM_NULLPTR );
#endif

  /* Insert the values */
  internal::insert( ext_vertex, n_IDs, values );
  internal::insert( ext_hex, n_IDs, values );

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( ext_vertex.getNumberOfValues(), n_IDs * vertex_stride);

  EXPECT_FALSE( ext_hex.empty() );
  EXPECT_EQ( ext_hex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( ext_hex.getNumberOfValues(), n_IDs * hex_stride);

  /* Check that the external ConnectivityArrays cannot insert any more. */
  EXPECT_DEATH_IF_SUPPORTED( ext_vertex.insert( values, 0 ), IGNORE_OUTPUT ); 
  EXPECT_DEATH_IF_SUPPORTED( ext_hex.insert( values, 0 ), IGNORE_OUTPUT ); 

  delete[] values;
  delete[] ext_vertex_values;
  delete[] ext_hex_values;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, NoIndirectionSidreInsert )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType n_IDs = 500 * 6;
  constexpr IndexType n_values = hex_stride * n_IDs;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ n_values ];
  for ( IndexType i = 0; i < n_values; ++i )
  {
    values[ i ] = i;
  }
  
  /* Create the sidre storage ConnectivityArrays to be tested. */
  ConnectivityArray< vertex_stride, vertex > 
                          sidre_vertex( root->createGroup( "vertex" ), "test" );
  EXPECT_TRUE( sidre_vertex.isInSidre() );
  EXPECT_FALSE( sidre_vertex.isExternal() );
  EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
  EXPECT_TRUE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< hex_stride, hex >
                          sidre_hex( root->createGroup( "hex" ), "test" );
  EXPECT_TRUE( sidre_hex.isInSidre() );
  EXPECT_FALSE( sidre_hex.isExternal() );
  EXPECT_EQ( sidre_hex.getGroup(), root->getGroup( "hex") );
  EXPECT_TRUE( sidre_hex.empty() );
  EXPECT_EQ( sidre_hex.getNumberOfIDs(), 0);

  /* Insert the values */
  internal::insert( sidre_vertex, n_IDs, values );
  internal::insert( sidre_hex, n_IDs, values );

  /* Check that the number of IDs is correct. */
  EXPECT_FALSE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( sidre_vertex.getNumberOfValues(), n_IDs * vertex_stride);

  EXPECT_FALSE( sidre_hex.empty() );
  EXPECT_EQ( sidre_hex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( sidre_hex.getNumberOfValues(), n_IDs * hex_stride);

  /* Check that restoring from sidre functions properly */
  sidre::Group* vertex_group = 
                        const_cast< sidre::Group* >( sidre_vertex.getGroup() );
  ConnectivityArray< vertex_stride, vertex > sidre_vertex_cpy( vertex_group );
  internal::check_equality( sidre_vertex, sidre_vertex_cpy, __LINE__ );

  sidre::Group* hex_group = 
                        const_cast< sidre::Group* >( sidre_hex.getGroup() );
  ConnectivityArray< hex_stride, hex > sidre_hex_cpy( hex_group );
  internal::check_equality( sidre_hex, sidre_hex_cpy, __LINE__ );

  delete[] values;

#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, IndirectionNativeInsert )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType n_IDs = 500 * 6;
  constexpr IndexType n_values = (vertex_stride + hex_stride) * n_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ n_values ];
  for ( IndexType i = 0; i < n_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ n_IDs + 1 ];
  CellTypes* types = new CellTypes[ n_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < n_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }
  
  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray< NEED_INDIRECTION, vertex > native_vertex;
  EXPECT_FALSE( native_vertex.isExternal() );
  EXPECT_FALSE( native_vertex.isInSidre() );
  EXPECT_TRUE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< NEED_INDIRECTION, mixed > native_mixed;
  EXPECT_FALSE( native_mixed.isExternal() );
  EXPECT_FALSE( native_mixed.isInSidre() );
  EXPECT_TRUE( native_mixed.empty() );
  EXPECT_EQ( native_mixed.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( native_mixed.getGroup(), AXOM_NULLPTR );
#endif

  /* Insert the values. */
  internal::insert( native_vertex, n_IDs, values, offsets );
  internal::insert( native_mixed, n_IDs, values, offsets, types );

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( native_vertex.empty() );
  EXPECT_EQ( native_vertex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( native_vertex.getNumberOfValues(), n_values);

  EXPECT_FALSE( native_mixed.empty() );
  EXPECT_EQ( native_mixed.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( native_mixed.getNumberOfValues(), n_values);

  delete[] values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, IndirectionExternalInsert )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType n_IDs = 500 * 6;
  constexpr IndexType n_values = (vertex_stride + hex_stride) * n_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ n_values ];
  for ( IndexType i = 0; i < n_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ n_IDs + 1 ];
  CellTypes* types = new CellTypes[ n_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < n_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }
  
  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[ n_values ];
  IndexType* ext_vertex_offsets = new IndexType[ n_IDs + 1 ];
  ConnectivityArray< NEED_INDIRECTION, vertex > 
                          ext_vertex( 0, ext_vertex_values, ext_vertex_offsets,
                                      n_IDs, n_values );
  EXPECT_TRUE( ext_vertex.isExternal() );
  EXPECT_FALSE( ext_vertex.isInSidre() );
  EXPECT_TRUE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), 0);

  IndexType* ext_mixed_values = new IndexType[ n_IDs * hex_stride ];
  IndexType* ext_mixed_offsets = new IndexType[ n_IDs + 1 ];
  CellTypes* ext_mixed_types = new CellTypes[ n_IDs ];
  ConnectivityArray< NEED_INDIRECTION, mixed > 
                          ext_mixed( 0, ext_mixed_values, ext_mixed_offsets, 
                                     ext_mixed_types, n_IDs, n_values );
  EXPECT_TRUE( ext_mixed.isExternal() );
  EXPECT_FALSE( ext_mixed.isInSidre() );
  EXPECT_TRUE( ext_mixed.empty() );
  EXPECT_EQ( ext_mixed.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
  EXPECT_EQ( ext_vertex.getGroup(), AXOM_NULLPTR );
  EXPECT_EQ( ext_mixed.getGroup(), AXOM_NULLPTR );
#endif

  /* Insert the values. */
  internal::insert( ext_vertex, n_IDs, values, offsets );
  internal::insert( ext_mixed, n_IDs, values, offsets, types );

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( ext_vertex.empty() );
  EXPECT_EQ( ext_vertex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( ext_vertex.getNumberOfValues(), n_values);

  EXPECT_FALSE( ext_mixed.empty() );
  EXPECT_EQ( ext_mixed.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( ext_mixed.getNumberOfValues(), n_values);

  /* Check that the external ConnectivityArrays cannot append any more. */
  EXPECT_DEATH_IF_SUPPORTED( ext_vertex.insert( values, 0, 1 ), IGNORE_OUTPUT ); 
  EXPECT_DEATH_IF_SUPPORTED( ext_mixed.insert( values,  0, 1, vertex ), 
                             IGNORE_OUTPUT ); 

  delete[] values;
  delete[] offsets;
  delete[] types;
  delete[] ext_vertex_values;
  delete[] ext_vertex_offsets;
  delete[] ext_mixed_values;
  delete[] ext_mixed_offsets;
  delete[] ext_mixed_types;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array, IndirectionSidreInsert )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType n_IDs = 500 * 6;
  constexpr IndexType n_values = (vertex_stride + hex_stride) * n_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ n_values ];
  for ( IndexType i = 0; i < n_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ n_IDs + 1 ];
  CellTypes* types = new CellTypes[ n_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < n_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  /* Create the sidre storage ConnectivityArrays to be tested. */
  ConnectivityArray< NEED_INDIRECTION, vertex > 
                          sidre_vertex( root->createGroup( "vertex" ), "test" );
  EXPECT_TRUE( sidre_vertex.isInSidre() );
  EXPECT_FALSE( sidre_vertex.isExternal() );
  EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
  EXPECT_TRUE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

  ConnectivityArray< NEED_INDIRECTION, mixed >
                          sidre_mixed( root->createGroup( "mixed" ), "test" );
  EXPECT_TRUE( sidre_mixed.isInSidre() );
  EXPECT_FALSE( sidre_mixed.isExternal() );
  EXPECT_EQ( sidre_mixed.getGroup(), root->getGroup( "mixed" ) );
  EXPECT_TRUE( sidre_mixed.empty() );
  EXPECT_EQ( sidre_mixed.getNumberOfIDs(), 0);

  /* Insert the values */
  internal::insert( sidre_vertex, n_IDs, values, offsets );
  internal::insert( sidre_mixed, n_IDs, values, offsets, types );

  /* Check that the number of IDs and values is correct. */
  EXPECT_FALSE( sidre_vertex.empty() );
  EXPECT_EQ( sidre_vertex.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( sidre_vertex.getNumberOfValues(), n_values);

  EXPECT_FALSE( sidre_mixed.empty() );
  EXPECT_EQ( sidre_mixed.getNumberOfIDs(), n_IDs);
  EXPECT_EQ( sidre_mixed.getNumberOfValues(), n_values);

  /* Check that restoring from sidre functions properly */
  sidre::Group* vertex_group = 
                        const_cast< sidre::Group* >( sidre_vertex.getGroup() );
  ConnectivityArray< NEED_INDIRECTION, vertex > sidre_vertex_cpy( vertex_group );
  internal::check_equality( sidre_vertex, sidre_vertex_cpy, __LINE__ );

  sidre::Group* mixed_group = 
                        const_cast< sidre::Group* >( sidre_mixed.getGroup() );
  ConnectivityArray< NEED_INDIRECTION, mixed > sidre_mixed_cpy( mixed_group );
  internal::check_equality( sidre_mixed, sidre_mixed_cpy, __LINE__ );

  delete[] values;
  delete[] offsets;
  delete[] types;

#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}


/*******************************************************************************
 *                          Capacity tests                                     *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, NoIndirectionNativeCapacity )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType n_IDs = 100;
  constexpr IndexType max_values = hex_stride * n_IDs;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }
  
  for ( double ratio = 0.75; ratio <= 3.0; ratio += 0.25 )
  {
    /* Create the native storage ConnectivityArrays to be tested. */
    ConnectivityArray< vertex_stride, vertex > native_vertex;
    native_vertex.setResizeRatio( ratio );
    EXPECT_FALSE( native_vertex.isExternal() );
    EXPECT_FALSE( native_vertex.isInSidre() );
    EXPECT_TRUE( native_vertex.empty() );
    EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

    ConnectivityArray< hex_stride, CellTypes::HEX > native_hex;
    native_hex.setResizeRatio( ratio );
    EXPECT_FALSE( native_hex.isExternal() );
    EXPECT_FALSE( native_hex.isInSidre() );
    EXPECT_TRUE( native_hex.empty() );
    EXPECT_EQ( native_hex.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
    EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
    EXPECT_EQ( native_hex.getGroup(), AXOM_NULLPTR );
#endif

    /* Append the values */
    if ( ratio < 1.0 )
    {
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( native_vertex, n_IDs, values ), 
        IGNORE_OUTPUT );
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( native_hex, n_IDs, values ), 
        IGNORE_OUTPUT );
    }
    else
    {
      internal::append_capacity( native_vertex, n_IDs, values );
      internal::append_capacity( native_hex, n_IDs, values );
    }
  }

  delete[] values;
}


//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, NoIndirectionSidreCapacity )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr IndexType n_IDs = 100;
  constexpr IndexType max_values = hex_stride * n_IDs;

  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }
  
  for ( double ratio = 0.75; ratio <= 3.0; ratio += 0.25 )
  {
    /* Create the sidre storage ConnectivityArrays to be tested. */
    ConnectivityArray< vertex_stride, vertex > 
                            sidre_vertex( root->createGroup( "vertex" ), "test" );
    sidre_vertex.setResizeRatio( ratio );
    EXPECT_TRUE( sidre_vertex.isInSidre() );
    EXPECT_FALSE( sidre_vertex.isExternal() );
    EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
    EXPECT_TRUE( sidre_vertex.empty() );
    EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

    ConnectivityArray< hex_stride, hex >
                            sidre_hex( root->createGroup( "hex" ), "test" );
    sidre_hex.setResizeRatio( ratio );
    EXPECT_TRUE( sidre_hex.isInSidre() );
    EXPECT_FALSE( sidre_hex.isExternal() );
    EXPECT_EQ( sidre_hex.getGroup(), root->getGroup( "hex") );
    EXPECT_TRUE( sidre_hex.empty() );
    EXPECT_EQ( sidre_hex.getNumberOfIDs(), 0);

    /* Append the values */
    if ( ratio < 1.0 )
    {
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( sidre_vertex, n_IDs, values ), 
        IGNORE_OUTPUT );
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( sidre_hex, n_IDs, values ), 
        IGNORE_OUTPUT );
    }
    else
    {
      internal::append_capacity( sidre_vertex, n_IDs, values );
      internal::append_capacity( sidre_hex, n_IDs, values );
    }

    root->destroyGroups();
  }

  delete[] values;
#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, IndirectionNativeCapacity )
{
  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType n_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * n_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ n_IDs + 1 ];
  CellTypes* types = new CellTypes[ n_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < n_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }
  
  for ( double ratio = 0.75; ratio <= 3.0; ratio += 0.25 )
  {
    /* Create the native storage ConnectivityArrays to be tested. */
    ConnectivityArray< NEED_INDIRECTION, vertex > native_vertex;
    native_vertex.setResizeRatio( ratio );
    EXPECT_FALSE( native_vertex.isExternal() );
    EXPECT_FALSE( native_vertex.isInSidre() );
    EXPECT_TRUE( native_vertex.empty() );
    EXPECT_EQ( native_vertex.getNumberOfIDs(), 0);

    ConnectivityArray< NEED_INDIRECTION, mixed > native_mixed;
    native_mixed.setResizeRatio( ratio );
    EXPECT_FALSE( native_mixed.isExternal() );
    EXPECT_FALSE( native_mixed.isInSidre() );
    EXPECT_TRUE( native_mixed.empty() );
    EXPECT_EQ( native_mixed.getNumberOfIDs(), 0);

#ifdef MINT_USE_SIDRE
    EXPECT_EQ( native_vertex.getGroup(), AXOM_NULLPTR );
    EXPECT_EQ( native_mixed.getGroup(), AXOM_NULLPTR );
#endif

    /* Append the values */
    if ( ratio < 1.0 )
    {
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( native_vertex, n_IDs, values, offsets ), 
        IGNORE_OUTPUT );
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( native_mixed, n_IDs, values, offsets, types ), 
        IGNORE_OUTPUT );
    }
    else
    {
      internal::append_capacity( native_vertex, n_IDs, values, offsets );
      internal::append_capacity( native_mixed, n_IDs, values, offsets, types );
    }
  }

  delete[] values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST( mint_connectivity_array_DeathTest, IndirectionSidreCapacity )
{
#ifdef MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr CellTypes vertex = CellTypes::VERTEX;
  constexpr IndexType vertex_stride = cell_info< vertex >::num_nodes;
  
  constexpr CellTypes hex = CellTypes::HEX;
  constexpr IndexType hex_stride = cell_info< hex >::num_nodes;

  constexpr CellTypes mixed = CellTypes::MIXED;

  constexpr IndexType n_IDs = 100;
  constexpr IndexType max_values = (vertex_stride + hex_stride) * n_IDs / 2;
  
  /* Allocate and populate the values buffer. */
  IndexType* values = new IndexType[ max_values ];
  for ( IndexType i = 0; i < max_values; ++i )
  {
    values[ i ] = i;
  }

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets = new IndexType[ n_IDs + 1 ];
  CellTypes* types = new CellTypes[ n_IDs ];
  offsets[ 0 ] = 0;
  for ( IndexType i = 0; i < n_IDs; ++i ) 
  {
    if ( i % 2 == 0 )
    {
      types[ i ] = vertex;
      offsets[ i + 1 ] = offsets[ i ] + vertex_stride;
    }
    else
    {
      types[ i ] = hex;
      offsets[ i + 1 ] = offsets[ i ] + hex_stride;
    }
  }

  for ( double ratio = 0.75; ratio <= 3.0; ratio += 0.25 )
  {
    /* Create the sidre storage ConnectivityArrays to be tested. */
    ConnectivityArray< NEED_INDIRECTION, vertex > 
                            sidre_vertex( root->createGroup( "vertex" ), "test" );
    sidre_vertex.setResizeRatio( ratio );
    EXPECT_TRUE( sidre_vertex.isInSidre() );
    EXPECT_FALSE( sidre_vertex.isExternal() );
    EXPECT_EQ( sidre_vertex.getGroup(), root->getGroup( "vertex" ) );
    EXPECT_TRUE( sidre_vertex.empty() );
    EXPECT_EQ( sidre_vertex.getNumberOfIDs(), 0);

    ConnectivityArray< NEED_INDIRECTION, mixed >
                            sidre_mixed( root->createGroup( "mixed" ), "test" );
    sidre_mixed.setResizeRatio( ratio );
    EXPECT_TRUE( sidre_mixed.isInSidre() );
    EXPECT_FALSE( sidre_mixed.isExternal() );
    EXPECT_EQ( sidre_mixed.getGroup(), root->getGroup( "mixed" ) );
    EXPECT_TRUE( sidre_mixed.empty() );
    EXPECT_EQ( sidre_mixed.getNumberOfIDs(), 0);

    /* Append the values */
    if ( ratio < 1.0 )
    {
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( sidre_vertex, n_IDs, values, offsets ), 
        IGNORE_OUTPUT );
      EXPECT_DEATH_IF_SUPPORTED( 
        internal::append_capacity( sidre_mixed, n_IDs, values, offsets, types ), 
        IGNORE_OUTPUT );
    }
    else
    {
      internal::append_capacity( sidre_vertex, n_IDs, values, offsets );
      internal::append_capacity( sidre_mixed, n_IDs, values, offsets, types );
    }

    root->destroyGroups();
  }

  delete[] values;
  delete[] offsets;
  delete[] types;

#else
  EXPECT_TRUE( true );
#endif  /* MINT_USE_SIDRE */
}


}   /* end namespace mint */
}   /* end namespace axom */


//------------------------------------------------------------------------------
using axom::slic::UnitTestLogger;

int main( int argc, char* argv[] )
{
  int result = 0;
  ::testing::InitGoogleTest( &argc, argv );
  UnitTestLogger logger;
  result = RUN_ALL_TESTS();
  return result;
}
