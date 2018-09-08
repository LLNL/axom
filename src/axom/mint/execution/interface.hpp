/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef MINT_EXECUTION_INTERFACE_HPP_
#define MINT_EXECUTION_INTERFACE_HPP_

#include "axom/config.hpp"      // compile-time definitions
#include "axom/core/Macros.hpp" // for AXOM_STATIC_ASSERT

#include "axom/mint/execution/xargs.hpp"                   // for xargs
#include "axom/mint/execution/internal/for_all_cells.hpp"  // for_all_cells()
#include "axom/mint/execution/internal/for_all_nodes.hpp"  // for_all_nodes()
#include "axom/mint/execution/policy.hpp"                  // for policy_traits

#include "axom/mint/mesh/Mesh.hpp"                         // for  mint::Mesh

// Slic includes
#include "axom/slic/interface/slic.hpp"                    // for SLIC macros

#ifdef AXOM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

// C/C++ includes
#include <type_traits> // for std::is_same

/*!
 * \file
 *
 * \brief Defines an execution interface for generic loop traversals and various
 *  mesh traversals, e.g., loop over nodes, cells, etc. The traversals may be
 *  executed in parallel, on the CPU, GPU, or other accelerator device, or
 *  in serial on the CPU according to the execution policy, which is supplied as
 *  a template argument to each of the traversal functions.
 *
 *  The general interface for the mesh traversal functions takes in a mint Mesh
 *  object and the loop body encapsulated in a lambda expression, conforming to
 *  the following template:
 *  \code
 *     for_all_[entity]< [exec_policy], [xargs] >( meshPtr, AXOM_LAMBDA(...) {
 *          // loop body
 *     } );
 *  \endcode
 *
 *  where:
 *
 *   * <b> [entity] </b> <br/>
 *     The entity suffix specifies the mesh entity being traversed, e.g., nodes,
 *     cells, faces, etc,.
 *
 *   * <b> [exec_policy] </b> <br />
 *     The execution policy indicates <em> how </em> and <em> where </em>
 *     the data corresponding to the specified mesh entity will be traversed.
 *     For example, an execution policy may indicate that the mesh traversal
 *     will be executed in parallel on the GPU, or CPU. A list of the
 *     currently supported execution policies and a brief description is given
 *     in policy.hpp
 *
 *   * <b> [AXOM_LAMBDA()] </b> <br />
 *     The AXOM_LAMBDA argument encapsulates the loop body, i.e., the kernel
 *     that is to be executed at each mesh entity. By default, the lambda
 *     expression takes the ID of the corresponding mesh entity as an argument.
 *     Additional lambda arguments may be specified by supplying the
 *     <em> xargs </em> template argument to the mesh traversal function.
 *
 *   * <b> [xargs] </b> <br />
 *     The xargs is an optional template argument to the mesh traversal function
 *     that specifies any additional arguments that the supplied lambda
 *     expression takes.A list of the available xargs types and a brief
 *     description is given in xargs.hpp
 *
 *  \note For parallel execution, the implementation relies heavily on the
 *   <a href="https://github.com/LLNL/RAJA"> RAJA </a> programming model
 *   abstraction layer.
 *
 *  \see policy.hpp
 *  \see xargs.hpp
 */

namespace axom
{
namespace mint
{

/// \name Generic Loop Traversal Functions
/// @{

/*!
 * \brief Iterates over the specified interval, \f$\mathcal{I}\f$:[begin,end-1],
 *  and executes the supplied kernel at each iteration.
 *
 * \param [in] begin the starting index of the iteration.
 * \param [in] end the total number of iterations
 * \param [in] kernel lambda expression consisting of the loop body
 *
 * \tparam ExecPolicy the execution policy
 * \tparam KernelType the type of the supplied lambda expression.
 *
 * \note The supplied lambda expression is expected to take the index of the
 *  iteration as an argument.
 *
 * Usage Example:
 * \code
 *
 *   constexpr mint::IndexType N = 10;
 *   double a[ N ];
 *   double b[ N ];
 *   double c[ N];
 *
 *   // other code
 *   ...
 *
 *   using exec = policy::parallel_cpu;
 *   mint::for_all< exec >( 0, N, AXOM_LAMBDA(mint::IndexType idx) {
 *     c[ idx ] = a[ idx ] + b[ idx ];
 *   } );
 *
 * \endcode
 *
 * \see policy.hpp
 */
template < typename ExecPolicy = policy::serial, typename KernelType >
inline void for_all( const mint::IndexType& begin,
                     const mint::IndexType& end,
                     KernelType&& kernel )
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT( policy_traits< ExecPolicy >::valid() );

#ifdef AXOM_USE_RAJA
  using exec_pol = typename policy_traits< ExecPolicy >::raja_exec_policy;
  RAJA::forall< exec_pol >( RAJA::RangeSegment(begin,end), kernel );
#else

  constexpr bool is_serial = std::is_same< ExecPolicy, policy::serial >::value;
  AXOM_STATIC_ASSERT( is_serial );
  for ( mint::IndexType i=begin ; i < end ; ++i )
  {
    kernel( i );
  }

#endif

}

/*!
 * \brief Iterates over the interval, \f$\mathcal{I}\f$:[0, N-1], and
 *  executes the supplied kernel at each iteration.
 *
 * \param [in] N the total number of iterations
 * \param [in] kernel lambda expression consisting of the loop body
 *
 * \tparam ExecPolicy the execution policy
 * \tparam KernelType the type of the supplied lambda expression.
 *
 * \note The supplied lambda expression is expected to take the index of the
 *  iteration as an argument.
 *
 * Usage Example:
 * \code
 *
 *   constexpr int N = 10;
 *   double a[ N ];
 *   double b[ N ];
 *   double c[ N];
 *
 *   // other code
 *   ...
 *
 *   using exec = policy::parallel_cpu;
 *   mint::for_all< exec >( N, AXOM_LAMBDA(mint::IndexType idx) {
 *     c[ idx ] = a[ idx ] + b[ idx ];
 *   } );
 *
 * \endcode
 *
 * \see policy.hpp
 */
template < typename ExecPolicy=policy::serial, typename KernelType >
inline void for_all( const mint::IndexType& N, KernelType&& kernel )
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT( policy_traits< ExecPolicy >::valid() );
  for_all< ExecPolicy >( 0, N, std::forward< KernelType >( kernel ) );
}

/// @}

/// \name Mesh Node Traversal Functions
/// @{

/*!
 * \brief Loops over the nodes of the given mesh.
 *
 * \param [in] m pointer to the mesh object
 * \param [in] kernel user-supplied lambda consisting of the kernel
 *
 * \pre m != nullptr
 *
 * \tparam ExecPolicy the execution policy, e.g., serial or parallel
 * \tparam ArgType object indicating the arguments to the kernel
 * \tparam KernelType
 *
 * Usage Example:
 * \code
 *
 *   for_all_nodes< exec >( m, AXOM_LAMDA( IndexType nodeIdx ) {
 *      foo[ nodeIdx ] = val;
 *      ...
 *   } );
 *
 *
 *   for_all_nodes< exec, xargs::ij >( m,
 *          AXOM_LAMDA( IndexType nodeIdx, IndexType i, IndexType j) {
 *     ...
 *   } );
 *
 *   for_all_nodes< exec, xargs::ijk >( m,
 *      AXOM_LAMDA( IndexType nodeIdx, IndexType i, IndexType j, IndexType k) {
 *     ...
 *   } );
 *
 *   for_all_nodes< exec, xargs::xy >( m,
 *           AXOM_LAMBA( IndexType nodeIdx, double x, double y) {
 *     ...
 *   } );
 *
 *   for_all_nodes< exec, xargs::xyz >( m,
 *           AXOM_LAMBA( IndexType nodeIdx, double x, double y, double z) {
 *     ...
 *   } );
 * \endcode
 *
 * \see policy.hpp
 * \see xargs.hpp
 */
template < typename ExecPolicy = policy::serial,
           typename ArgType    = xargs::index,
           typename KernelType >
inline void for_all_nodes( const mint::Mesh* m, KernelType&& kernel )
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT( policy_traits< ExecPolicy >::valid() );
  AXOM_STATIC_ASSERT( xargs_traits< ArgType >::valid ());

  // run-time sanity checks
  SLIC_ASSERT( m != nullptr );

  // dispatch
  internal::for_all_nodes< ExecPolicy >(
    ArgType(), m, std::forward< KernelType >( kernel ) );
}

/// @}

/// \name Mesh Cell Traversal Functions
/// @{

/*!
 * \brief Loops over all the cells of a given mesh.
 *
 * \param [in] m pointer to the mesh object.
 * \param [in] kernel user-supplied kernel to execute on each cell.
 *
 * \pre m != nullptr
 *
 * \tparam ExecPolicy the execution policy, e.g., serial or parallel
 * \tparam ArgType object indicating the arguments to the kernel
 *
 * Usage Example:
 * \code
 *
 *   for_all_cells< exec >( m, AXOM_LAMDA( IndexType cellIdx ) {
 *      foo[ cellIdx ] = val;
 *      ...
 *   } );
 *
 *
 *   for_all_cells< exec, xargs::ij >( m,
 *          AXOM_LAMDA( IndexType cellIdx, IndexType i, IndexType j) {
 *     ...
 *   } );
 *
 *   for_all_cells< exec, xargs::ijk >( m,
 *      AXOM_LAMDA( IndexType cellIdx, IndexType i, IndexType j, IndexType k) {
 *     ...
 *   } );
 *
 *  for_all_cells< exec, xargs::nodeids >( m,
 *     AXOM_LAMBDA( IndexType cellIdx, const IndexType* nodeIds, IndexType N )
 *  } );
 *
 * \endcode
 *
 * \see policy.hpp
 * \see xargs.hpp
 */
template < typename ExecPolicy = policy::serial,
           typename ArgType    = xargs::index,
           typename KernelType >
inline void for_all_cells( const mint::Mesh* m, KernelType&& kernel )
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT( policy_traits< ExecPolicy >::valid() );
  AXOM_STATIC_ASSERT( xargs_traits< ArgType >::valid() );

  // run-time sanity checks
  SLIC_ASSERT( m != nullptr );

  // dispatch
  internal::for_all_cells< ExecPolicy >(
    ArgType(), m, std::forward< KernelType >( kernel )  );
}

/// @}

/// \name Synchronization
/// @{

/*!
 * \brief Synchronizes all execution threads of the specified execution policy.
 *
 * \tparam ExecPolicy the mint execution policy
 */
/// @{
template < typename ExecPolicy >
inline void synchronize( )
{
  // compile-time sanity checks
  AXOM_STATIC_ASSERT( policy_traits< ExecPolicy >::valid() );

#ifdef AXOM_USE_RAJA
  using sync_pol = typename policy_traits< ExecPolicy >::raja_sync_policy;
  RAJA::synchronize< sync_pol >( );
#endif
}

template < >
inline void synchronize< policy::serial >( ) {
  return;
}
/// @}

/// @}

}
}
#endif /* MINT_EXECUTION_INTERFACE_HPP_ */
