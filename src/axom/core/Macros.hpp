// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file AxomMacros.hpp
 *
 * \brief Contains several useful macros for the axom project
 *
 */

#ifndef AXOM_MACROS_HPP_
#define AXOM_MACROS_HPP_

#include "axom/config.hpp"
#include <cassert>  // for assert()

/*!
 * \def AXOM_DEVICE
 * \def AXOM_HOST_DEVICE
 *
 * \brief CUDA host/device macros for decorating functions/lambdas
 *
 * \note These will expand to the corresponding CUDA decorations when
 *  compiled with -DAXOM_USE_CUDA
 */
#if defined(__CUDACC__)
  #define AXOM_DEVICE __device__
  #define AXOM_HOST_DEVICE __host__ __device__
  #define AXOM_HOST __host__
#else
  #define AXOM_DEVICE
  #define AXOM_HOST_DEVICE
  #define AXOM_HOST
#endif

/*
 * \def AXOM_STRINGIFY
 *
 * \brief Helper Macro to specify strings inside pragmas.
 */
#define AXOM_STRINGIFY(x) AXOM_DO_STRINGIFY(x)
#define AXOM_DO_STRINGIFY(x) #x

/*
 * \def AXOM_PRAGMA
 *
 * \brief Macro used to specify pragma directive
 */
#ifdef _WIN32
  #define AXOM_PRAGMA(x) __pragma(x)
#else
  #define AXOM_PRAGMA(x) _Pragma(AXOM_STRINGIFY(x))
#endif

/*
 * \def AXOM_SUPPRESS_HD_WARN
 *
 * \brief Macro used to silence __host__ __device__ compiler warnings
 *  when calling __host__ function from __host__ __device__ function.
 */
#if defined(__CUDACC__)
  #define AXOM_SUPPRESS_HD_WARN AXOM_PRAGMA(nv_exec_check_disable)
#else
  #define AXOM_SUPPRESS_HD_WARN
#endif

/*!
 * \def AXOM_LAMBDA
 *
 * \brief Convenience macro used for lambda capture by value.
 * \note When CUDA is used, the macro always expands to a host/device lambda.
 *
 * \warning When compiling with CUDA, host/device lambdas incur a significant
 *  penalty on the CPU code. The way NVCC implements host/device lambdas
 *  prevents the compiler from proper in-lining them. When CUDA is enabled use
 *  the parallel_gpu execution policy or opt to turn off CUDA if the application
 *  is making more use of the parallel_cpu and serial execution policies.
 */
#ifdef AXOM_USE_CUDA
  #define AXOM_LAMBDA [=] AXOM_HOST_DEVICE
  #define AXOM_DEVICE_LAMBDA [=] AXOM_DEVICE
  #define AXOM_HOST_LAMBDA [=] AXOM_HOST
#else
  #define AXOM_LAMBDA [=]
  #define AXOM_DEVICE_LAMBDA [=]
  #define AXOM_HOST_LAMBDA [=]
#endif

/*!
 * \def AXOM_CUDA_TEST
 *
 * \brief Convenience macro used for a gtest that uses cuda.
 */
#if defined(AXOM_USE_CUDA)
  #define AXOM_CUDA_TEST(X, Y)         \
    static void cuda_test_##X##Y();    \
    TEST(X, Y) { cuda_test_##X##Y(); } \
    static void cuda_test_##X##Y()
#else
  #define AXOM_CUDA_TEST(X, Y) TEST(X, Y)
#endif

/*!
 * \def AXOM_DEVICE_CODE
 *
 * \brief Convenience macro used for kernel code
 */
#if defined(__CUDA_ARCH__)
  #define AXOM_DEVICE_CODE
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
 * \brief This macro wraps the compile time static_assert functionality
 *  so you don't have to provide a message.
 */
#define AXOM_STATIC_ASSERT(cond) static_assert(cond, #cond)
#define AXOM_STATIC_ASSERT_MSG(cond, MSG) static_assert(cond, MSG)

/*!
 *
 * \def AXOM_UNUSED_VAR(x)
 * \brief Macro used to silence compiler warnings about variables
 *        that are defined but not used.
 * \note The intent is to use this macro for variables that are only used
 *       within compiler defines (e.g. in debug assertions). For example:
 * \code
 *
 *  double myVar = ...
 *  AXOM_UNUSED_VAR(myVar);
 *
 *  // code emits the following warning in release builds
 *  // if extra warnings are enabled and this macro is not called
 *  // warning: unused variable 'myVar' [-Wunused-variable]
 *
 *  SLIC_ASSERT(myVar > 0)
 *
 * \endcode
 */
#define AXOM_UNUSED_VAR(_x) static_cast<void>(_x)

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
  #define AXOM_DEBUG_PARAM(_x) _x
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
#define DISABLE_DEFAULT_CTOR(className) className() = delete

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
#define DISABLE_COPY_AND_ASSIGNMENT(className) \
  className(const className&) = delete;        \
  className& operator=(const className&) = delete

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
#define DISABLE_MOVE_AND_ASSIGNMENT(className) \
  className(const className&&) = delete;       \
  className& operator=(const className&&) = delete

#endif /* AXOM_MACROS_HPP_ */
