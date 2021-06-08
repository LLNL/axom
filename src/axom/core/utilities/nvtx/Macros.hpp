// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NVTX_MACROS_HPP_
#define AXOM_NVTX_MACROS_HPP_

#include "axom/core/utilities/nvtx/Range.hpp"

/*!
 * \file
 *
 * \brief Defines NVTX Macros that can be used to annotate functions and
 *  sections of the code. These macros use the NVIDIA Tools Extension library
 *  to provide additional information to NVIDIA performance tools, e.g.,
 *  nvprof, nvvp, Nsight, thereby, facilitate developers in performance
 *  evaluation.
 *
 */

/// \name AXOM NVTX Macros
///@{

/*!
 * \def AXOM_NVTX_SECTION
 *
 * \brief The AXOM_NVTX_SECTION macro is used to annotate sections of code
 *
 * \note In contrast to the AXOM_NVTX_FUNCTION macro, the AXOM_NVTX_SECTION
 *   macro is used to annotate sections of code, at a much finer granularity,
 *   within a given function.
 *
 * \warning Variables declared within a given AXOM_NVTX_SECTION are only defined
 *  within the scope of the AXOM_NVTX_SECTION.
 *
 * \warning An AXOM_NVTX_SECTION cannot be called in a nested fashion, i.e.,
 *  within another AXOM_NVTX_SECTION
 *
 * \note You may have multiple AXOM_NVTX_SECTION defined within a function and
 *  this macro can be used in conjunction with the AXOM_NVTX_FUNCTION macro.
 *
 * \Usage Example:
 * \code
 *
 *   void foo( )
 *   {
 *     AXOM_NVTX_FUNCTION( "foo"" );
 *
 *     // STEP 0: Run kernel A
 *     AXOM_NVTX_SECTION( "kernelA",
 *
 *        axom::for_all( 0, N, AXOM_LAMBDA(axom::IndexType i)
 *        {
 *          ..
 *        } );
 *
 *     ); // END NVTX SECTION for kernel A
 *
 *     // STEP 1: Run kernel B
 *     AXOM_NVTX_SECTION( "kernelB",
 *
 *        axom::for_all( 0, N, AXOM_LAMBDA(axom::IndexType i)
 *        {
 *          ...
 *        } );
 *
 *     ); // END NVTX SECTION for kernel B
 *
 *   }
 * \endcode
 *
 */
#if defined(AXOM_USE_ANNOTATIONS) && defined(AXOM_USE_CUDA)
  #define AXOM_NVTX_SECTION(__name__, ...) \
    do                                     \
    {                                      \
      axom::nvtx::Range r(__name__);       \
      __VA_ARGS__                          \
    } while(false)
#else
  #define AXOM_NVTX_SECTION(__name__, ...) \
    do                                     \
    {                                      \
      __VA_ARGS__                          \
    } while(false)
#endif

/*!
 * \def AXOM_NVTX_FUNCTION( name )
 *
 * \brief The AXOM_NVTX_FUNCTION macro is used to annotate a function.
 * \param [in] name a user-supplied name that will be given to the range.
 *
 * \note Typically, the AXOM_NVTX_FUNCTION macro is placed in the beginning of
 *  the function to annotate.
 *
 * \warning The AXOM_NVTX_FUNCTION can be called once within a (function) scope.
 *
 * Usage Example:
 * \code
 *   void foo( )
 *   {
 *     AXOM_NVTX_FUNCTION( "foo" );
 *     ...
 *   }
 * \endcode
 *
 */
#if defined(AXOM_USE_ANNOTATIONS) && defined(AXOM_USE_CUDA)
  #define AXOM_NVTX_FUNCTION(__name__) axom::nvtx::Range __func_range(__name__)
#else
  #define AXOM_NVTX_FUNCTION(__name__)
#endif

///@}

#endif /* AXOM_NVTX_MACROS_HPP_ */
