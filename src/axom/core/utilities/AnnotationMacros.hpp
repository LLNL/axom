// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_ANNOTATION_MACROS_HPP_
#define AXOM_ANNOTATION_MACROS_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_CALIPER
  #include "axom/core/utilities/nvtx/interface.hpp"
#endif

/*!
 * \def AXOM_PERF_MARK_FUNCTION( name )
 * 
 * \brief The AXOM_PERF_MARK_FUNCTION is used to annotate a function.
 * 
 * \param [in] name a user-supplied name to annotate the function.
 * 
 * \note Typically, the AXOM_PERF_MARK_FUNCTION is placed in the beginning of
 *  the function to annotate.
 * 
 * 
 * \warning The AXOM_PERF_MARK_FUNCTION can only be called once within a given
 *  (function) scope.
 * 
 * 
 * Usage Example:
 * \code
 *   void foo( )
 *   {
 *     AXOM_PERF_MARK_FUNCTION( "foo" );
 *     ... 
 *   }
 * \endcode
 */
#if defined(AXOM_USE_ANNOTATIONS) && defined(AXOM_USE_CALIPER)
  #error "Support for Caliper has not yet been implemented in Axom!"
#elif defined(AXOM_USE_ANNOTATIONS)
  #define AXOM_PERF_MARK_FUNCTION(__func_name__) \
    AXOM_NVTX_FUNCTION(__func_name__)
#else
  #define AXOM_PERF_MARK_FUNCTION(__func_name__)
#endif

/*!
 * \def AXOM_PERF_MARK_SECTION( name )
 *  
 * \brief The AXOM_PERF_MARK_SECTION macro is used to annotate sections of code.
 * 
 * \note In contrast to the AXOM_PERF_MARK_FUNCTION, the AXOM_PERF_MARK_SECTION
 *  macro is used to annotate sections of code at a much finer granularity
 *  within a given function and it may be use in conjunction with the 
 *  AXOM_PERF_MARK_FUNCTION macro.
 * 
 * \warning Variables declared within an AXOM_PERF_MARK_SECTION are only defined
 *  within the scope of the annotated section.
 * 
 * \warning The AXOM_PERF_MARK_SECTION macro may not be called in a nested
 *  fashion, i.e., within another AXOM_PERF_MARK_SECTION
 *
 * Usage Example:
 * \code
 *   void foo( )
 *   {
 *      AXOM_PERF_MARK_FUNCTION( "foo" );
 *      
 *      AXOM_PERF_MARK_SECTION( "kernel_A",
 *        axom::for_all( 0, N, AXOM_LAMBDA(axom::IndexType idx)
 *        { 
 *          ...
 *        } );
 *      );
 * 
 *      AXOM_PERF_MARK_SECTION( "kernel_B",
 *        axom::for_all( 0, N, AXOM_LAMBDA(axom::IndexType idx)
 *        {
 *          ...
 *        } );
 *      );
 * 
 *   }
 * \endcode
 * 
 */
#if defined(AXOM_USE_ANNOTATIONS) && defined(AXOM_USE_CALIPER)
  #error "Support for Caliper has not yet been implemented in Axom!"
#elif defined(AXOM_USE_ANNOTATIONS)
  #define AXOM_PERF_MARK_SECTION(__name__, ...) \
    AXOM_NVTX_SECTION(__name__, __VA_ARGS__)
#else
  #define AXOM_PERF_MARK_SECTION(__name__, ...) \
    do                                          \
    {                                           \
      __VA_ARGS__                               \
    } while(false)
#endif

#endif
