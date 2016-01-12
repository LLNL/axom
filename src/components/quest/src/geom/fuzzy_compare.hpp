/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \brief Templated routines for comparing values within some tolerance.
 *
 * \date Dec 17, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef FUZZY_COMPARE_HPP_
#define FUZZY_COMPARE_HPP_

// C/C++ includes
#include <cmath>
#include <limits>

namespace quest
{
namespace math
{

/*!
 *******************************************************************************
 * \brief Checks if two values A,B are the same within some tolerance.
 * \param [in] A user-supplied value to check
 * \param [in] B user-supplied value to check
 * \param [in] tol tolerance used. Default is set to machine epsilon.
 * \return true \iff A is equal to B, otherwise false.
 *******************************************************************************
 */
template < typename real >
bool fuzzy_compare( const real& A,
                    const real& B,
                    const real& tol=std::numeric_limits< real >::epsilon() )
{
   return( std::fabs(A-B) < tol );
}


} /* end math namespace */

} /* end quest namespace */

#endif /* FUZZY_COMPARE_HPP_ */
