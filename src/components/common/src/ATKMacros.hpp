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
 * \def DISABLE_COPY_AND_ASSIGNMENT(className)
 * \brief Macro to disable copy and assignment operations for the given class.
 *******************************************************************************
 */
#ifdef USE_CXX11
#define DISABLE_COPY_AND_ASSIGNMENT(className)                                \
  className( const className& ) = delete;                                     \
  className& operator=(const className&) = delete;
#else
#define DISABLE_COPY_AND_ASSIGNMENT(className)                                \
  className( const className& );                                              \
  className& operator=( const className& );
#endif


#endif /* ATKMACROS_HPP_ */
