// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 *
 * \brief Evaluates a Bezier Curve using deCasteljau's algorithm
 *
 */

#ifndef EVAL_BEZIER_HPP_ 
#define EVAL_BEZIER_HPP_

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/core/numerics/Matrix.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"

namespace axom
{
namespace primal
{

/*!
 * \brief Evaluates a Bezier Curve at particular parameter value between 0 and 1
 *
 * \param [in] bCurve BezierCurve to be evaluated
 * \param [in] t parameter value between 0 and 1 at which to evaluate
 * \return p the value of the Bezier Curve at t
 *
 */

template < typename T, int NDIMS >
Point< T, NDIMS > eval_bezier(const BezierCurve< T, NDIMS > bcurve, T t)
{
  Point< T, NDIMS > ptval(0.0,NDIMS);
  int ord = bcurve.getOrder();
  T dCarray[ord+1];
  
  for ( int i=0; i < NDIMS; i++)
  {
     for ( int j=0 ; j <= ord ; j++)
     {
       dCarray[j] = bcurve[j][i];
     }
     for ( int j=1 ; j <= ord ; j++)
     {
       for ( int k=0 ; k <= ord-j ; k++)
       {
         dCarray[k]=(1-t)*dCarray[k]+t*dCarray[k+1];  
       }        
     }
     ptval[i]=dCarray[0];
  }

  return ptval;
}

template < typename T, int NDIMS >
void split_bezier(const BezierCurve< T, NDIMS > bcurve, T t, BezierCurve< T, NDIMS >& c1, BezierCurve< T, NDIMS>& c2)
{ 
  int ord = bcurve.getOrder();
  //SLIC_ASSERT( (ord == c1.getOrder()) && (ord == c2.getOrder()) );
  T* dCarray = new T[NDIMS*2*(ord+1)];
  for ( int i=0; i < NDIMS; i++)
  {
     for ( int j=0 ; j <= ord ; j++)
     {
       dCarray[i*(2*(ord)+1)+j] = bcurve[j][i];
     }
     for ( int j=1 ; j <= ord ; j++)
     {
       dCarray[i*(2*(ord)+1)+(ord)+j]=dCarray[0];
       for ( int k=0 ; k <= ord-j ; k++)
       {
         dCarray[i*(2*(ord)+1)+k]=(1-t)*dCarray[i*(2*(ord)+1)+k]+t*dCarray[i*(2*(ord)+1)+k+1];
       }       
     }
  }
     c1=BezierCurve< T, NDIMS> (dCarray,ord+1);
     c2=BezierCurve< T, NDIMS> (dCarray+NDIMS*(ord),ord+1);
     delete [] dCarray;
     return;
}

} /* namespace primal */
} /* namespace axom */

#endif /* EVAL_BEZIER_HPP_ */
