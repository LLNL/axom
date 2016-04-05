/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file ATKMacros.hpp
 *
 * \date Jul 28, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef ATKMACROS_HPP_
#define ATKMACROS_HPP_

/*!
 *******************************************************************************
 * \def ATK_NOT_USED(x)
 * \brief Macro used to silence compiler warnings in methods with unused
 *  arguments.
 * \note The intent is to use this macro in the function signature. For example:
 * \code
 *
 *  void my_function(int x, int ATK_NOT_USED(y))
 *  {
 *    // my implementation
 *  }
 *
 * \endcode
 *******************************************************************************
 */
#define ATK_NOT_USED(x)




/*!
 *******************************************************************************
 * \def ATK_DEBUG_VARIABLE(x)
 * \brief Macro used to silence compiler warnings about variables
 *        that are defined but not used.
 * \note The intent is to use this macro for variables that are only used
 *       for debugging purposes. For example:
 * \code
 *
 *  double myVar = ...
 *  ATK_DEBUG_VARIABLE(myVar)   // code will emit the following warning if extra
 *                              // warnings are enabled and macro is not called
 *                              // warning: unused variable 'myVar' [-Wunused-variable]
 *  SLIC_ASSERT(myVar > 0)
 *
 * \endcode
 *******************************************************************************
 */
#define ATK_DEBUG_VARIABLE(x)   if (0 && &x == &x){}


/*!
 *******************************************************************************
 * \def ATK_DEBUG_PARAM(x)
 * \brief Macro used to silence compiler warnings about parameters
 *        that are used in debug code but not in release code.
 * \note Default values are ok
 * \code
 *
 *  void my_function(int x, int ATK_DEBUG_PARAM(y))
 *  {
 *    // my implementation
 *    SLIC_ASSERT(y > 0)
 *  }
 *
 * \endcode
 *******************************************************************************
 */
#ifdef ATK_DEBUG
 #define ATK_DEBUG_PARAM(_x)  _x
#else
 #define ATK_DEBUG_PARAM(_x)
#endif


/*!
 *******************************************************************************
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
 *******************************************************************************
 */
#ifdef USE_CXX11
#define DISABLE_COPY_AND_ASSIGNMENT(className)                                \
  className( const className & ) = delete;                                     \
  className& operator=(const className&) = delete
#else
#define DISABLE_COPY_AND_ASSIGNMENT(className)                                \
  className( const className & );                                              \
  className& operator=( const className& )
#endif


#endif /* ATKMACROS_HPP_ */
