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
#include "axom/config.hpp"           // defines AXOM_USE_CXX11

/*!
 *
 * \file AxomMacros.hpp
 *
 * \brief Contains several useful macros for the axom project
 *
 */

#ifndef AXOM_MACROS_HPP_
#define AXOM_MACROS_HPP_

/*!
 * \def AXOM_DEVICE
 * \def AXOM_HOST_DEVICE
 *
 * \brief CUDA host/device macros for decorating functions/lambdas
 *
 * \note These will expand to the corresponding CUDA decorations when
 *  compiled with -DAXOM_USE_CUDA
 */
#ifdef AXOM_USE_CUDA
#define AXOM_DEVICE __device__
#define AXOM_HOST_DEVICE __host__ __device__
#else
#define AXOM_DEVICE
#define AXOM_HOST_DEVICE
#endif

/*!
 * \def AXOM_LAMBDA
 *
 * \brief Convenience macro used for lambda capture by value.
 * \note When CUDA is used, the macro always expands to a device lambda.
 */
#ifdef AXOM_USE_CUDA
#define AXOM_LAMBDA [=] AXOM_DEVICE
#else
#define AXOM_LAMBDA [=]
#endif

/*!
 *
 * \def AXOM_NOT_USED(x)
 * \brief Macro used to silence compiler warnings in methods with unused
 *  arguments.
 * \note The intent is to use this macro in the function signature. For example:
 * \code
 *
 *  void my_function(int x, int AXOM_NOT_USED(y))
 *  {
 *    // my implementation
 *  }
 *
 * \endcode
 */
#define AXOM_NOT_USED(x)

/*!
 *
 * \def AXOM_STATIC_ASSERT(cond)
 * \def AXOM_STATIC_ASSERT_MSG(cond, MSG)
 *
 * \brief This macro wraps C++11 compile time static_assert functionality. Used
 *  for backwards compatibility with non C++11 compilers.
 */
#ifdef AXOM_USE_CXX11
#define AXOM_STATIC_ASSERT( cond ) static_assert( cond, #cond )
#define AXOM_STATIC_ASSERT_MSG( cond, MSG ) static_assert( cond, MSG )
#else
#define AXOM_STATIC_ASSERT( cond )
#define AXOM_STATIC_ASSERT_MSG( cond, MSG )
#endif

/*!
 *
 * \def AXOM_DEBUG_VAR(x)
 * \brief Macro used to silence compiler warnings about variables
 *        that are defined but not used.
 * \note The intent is to use this macro for variables that are only used
 *       for debugging purposes (e.g. in debug assertions). For example:
 * \code
 *
 *  double myVar = ...
 *  AXOM_DEBUG_VAR(myVar);
 *
 *  // code emits the following warning in release builds
 *  // if extra warnings are enabled and this macro is not called
 *  // warning: unused variable 'myVar' [-Wunused-variable]
 *
 *  SLIC_ASSERT(myVar > 0)
 *
 * \endcode
 */
#define AXOM_DEBUG_VAR(_x)   static_cast<void>(_x)

/*!
 * \def AXOM_DEBUG_PARAM(x)
 * \brief Macro used to silence compiler warnings about parameters
 *        that are used in debug code but not in release code.
 * \note Default values are ok
 * \code
 *
 *  void my_function(int x, int AXOM_DEBUG_PARAM(y))
 *  {
 *    // my implementation
 *    SLIC_ASSERT(y > 0)
 *  }
 *
 * \endcode
 */
#ifdef AXOM_DEBUG
 #define AXOM_DEBUG_PARAM(_x)  _x
#else
 #define AXOM_DEBUG_PARAM(_x)
#endif

/*!
 * \def DISABLE_DEFAULT_CTOR(className)
 * \brief Macro to disable default constructor for the given class.
 * \note This macro should only be used within the private section of a class,
 *  as indicated in the example below.
 *
 * \code
 *
 *   class Foo
 *   {
 *   public:
 *
 *       // Public methods here
 *
 *   private:
 *      DISABLE_DEFAULT_CTOR(Foo);
 *   };
 *
 * \endcode
 */
#ifdef AXOM_USE_CXX11
#define DISABLE_DEFAULT_CTOR(className)                      \
  className( ) = delete;
#else
#define DISABLE_DEFAULT_CTOR(className)                      \
  className( );
#endif

/*!
 * \def DISABLE_COPY_AND_ASSIGNMENT(className)
 * \brief Macro to disable copy and assignment operations for the given class.
 * \note This macro should only be used within the private section of a class,
 *  as indicated in the example below.
 *
 * \code
 *
 *   class Foo
 *   {
 *   public:
 *      Foo();
 *      ~Foo();
 *
 *       // Other methods here
 *
 *   private:
 *      DISABLE_COPY_AND_ASSIGNMENT(Foo);
 *   };
 *
 * \endcode
 */
#ifdef AXOM_USE_CXX11
#define DISABLE_COPY_AND_ASSIGNMENT(className)                                 \
  className( const className & ) = delete;                                     \
  className& operator=(const className&) = delete
#else
#define DISABLE_COPY_AND_ASSIGNMENT(className)                                 \
  className( const className & );                                              \
  className& operator=( const className& )
#endif

/*!
 * \def DISABLE_MOVE_AND_ASSIGNMENT(className)
 * \brief Macro to disable move constructor and move assignment operations for
 * the given class.
 * \note This macro should only be used within the private section of a class,
 *  as indicated in the example below.
 *
 * \code
 *
 *   class Foo
 *   {
 *   public:
 *      Foo();
 *      ~Foo();
 *
 *       // Other methods here
 *
 *   private:
 *      DISABLE_MOVE_AND_ASSIGNMENT(Foo);
 *   };
 *
 * \endcode
 */
#ifdef AXOM_USE_CXX11
#define DISABLE_MOVE_AND_ASSIGNMENT(className)                                \
  className( const className&& ) = delete;                                    \
  className& operator=(const className&&) = delete
#else
#define DISABLE_MOVE_AND_ASSIGNMENT(className)
#endif


#endif /* AXOM_MACROS_HPP_ */
