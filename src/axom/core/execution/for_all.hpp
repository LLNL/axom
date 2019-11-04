// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_CORE_EXECUTION_FOR_ALL_HPP_
#define AXOM_CORE_EXECUTION_FOR_ALL_HPP_

#include "axom/config.hpp"                         /* compile time defs */
#include "axom/core/execution/execution_space.hpp" /* execution_space traits */
#include "axom/core/Macros.hpp"                    /* for axom Macros */
#include "axom/core/Types.hpp"                     /* for axom::IndexType */

// C/C++ includes
#include <type_traits> // for std::is_same()

namespace axom
{

/// \name Generic Loop Traversal Functions
/// @{

/*!
 * \brief
 *
 * \param begin
 * \param end
 * \param kernel
 */
template < typename ExecSpace, typename KernelType >
inline void for_all( const IndexType& begin, const IndexType& end,
                     KernelType&& kernel ) noexcept
{
  AXOM_STATIC_ASSERT( execution_space< ExecSpace >::valid() );

#ifdef AXOM_USE_RAJA

  using loop_exec = typename execution_space< ExecSpace >::loop_policy;
  RAJA::forall< loop_exec >( RAJA::RangeSegment( begin, end ),
                             std::forward< KernelType >( kernel ) );

#else

  constexpr bool is_serial = std::is_same< ExecSpace, SEQ_EXEC >::value;
  AXOM_STATIC_ASSERT( is_serial );
  for ( IndexType i=begin; i < end; ++i )
  {
    kernel( i );
  }

#endif

}

/*!
 * \brief
 *
 * \param N
 * \param kernel
 *
 */
template < typename ExecSpace, typename KernelType >
inline void for_all( const IndexType& N, KernelType&& kernel ) noexcept
{
  AXOM_STATIC_ASSERT( execution_space< ExecSpace >::valid() );
  for_all< ExecSpace >( 0, N, std::forward< KernelType >( kernel ) );
}

/// @}

} /* namespace axom */

#endif /* AXOM_CORE_EXECUTION_FOR_ALL_HPP_ */
