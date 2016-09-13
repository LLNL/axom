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
 * \file Triangle.hpp
 *
 * \date Dec 9, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef TRIANGLE_HPP_
#define TRIANGLE_HPP_

#include "quest/Point.hpp"
#include "quest/Vector.hpp"
#include "common/Utilities.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"

namespace quest
{

// Forward declare the templated classes and operator functions
template<typename T, int DIM> class Triangle;

/*!
 * \brief Overloaded output operator for triangles
 */
template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const Triangle<T,DIM> & tri);



template < typename T, int DIM >
class Triangle
{
public:
    typedef Point<T,DIM>  PointType;
    typedef Vector<T,DIM> VectorType;

    enum {
        NUM_TRI_VERTS = 3
    };

public:

  /*!
   *****************************************************************************
   * \brief Default constructor. Creates a degenerate triangle.
   *****************************************************************************
   */
  Triangle() { }

  /*!
   *****************************************************************************
   * \brief Custom Constructor. Creates a triangle from the 3 points A,B,C.
   * \param [in] A point instance corresponding to vertex A of the triangle.
   * \param [in] B point instance corresponding to vertex B of the triangle.
   * \param [in] C point instance corresponding to vertex C of the triangle.
   *****************************************************************************
   */
  Triangle( const PointType& A,
            const PointType& B,
            const PointType& C )
        : m_A (A), m_B(B), m_C(C)    {}


  /*!
   *****************************************************************************
   * \brief Destructor
   *****************************************************************************
   */
   ~Triangle() { }

  /*!
   *****************************************************************************
   * \brief Returns a reference to vertex A of the triangle.
   * \return A reference to vertex A of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& A( ) { return m_A; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex A of the triangle.
   * \return A const reference to vertex A of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& A() const { return m_A; };

  /*!
   *****************************************************************************
   * \brief Returns a reference to vertex B of the triangle.
   * \return B reference to vertex B of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& B( ) { return m_B; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex B of the triangle.
   * \return B const reference to vertex B of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& B() const { return m_B; };

  /*!
   *****************************************************************************
   * \brief Returns a reference to vertex C of the triangle.
   * \return C reference to vertex C of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& C( ) { return m_C; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex C of the triangle.
   * \return C const reference to vertex C of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& C() const { return m_C; };


  /**
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   */
  PointType& operator[](int idx)
  {
      SLIC_ASSERT(idx >=0 && idx < NUM_TRI_VERTS);

      return (idx == 0)
              ? m_A
              : (idx == 1)? m_B : m_C;
  }

  /**
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   */
  const PointType& operator[](int idx) const
  {
      SLIC_ASSERT(idx >=0 && idx < NUM_TRI_VERTS);

      return (idx == 0)
              ? m_A
              : (idx == 1)? m_B : m_C;
  }

  /*!
   * \brief Returns the normal of the triangle (not normalized)
   * \pre This function is only valid when DIM = 3
   * \return The normal vector to the triangle (when DIM==3), the zero vector otherwise
   */
  VectorType normal() const
  {
      SLIC_CHECK_MSG(DIM==3, "quest::Triangle::normal() is only valid when dimension is 3.");

      return (DIM==3)
              ? VectorType::cross_product( VectorType(m_A,m_B), VectorType(m_A,m_C))
              : VectorType();
  }

  /**
   * \brief Returns the area of the triangle
   * \pre Only defined when dimension DIM is 2 or 3
   */
  double area() const
  {
      SLIC_CHECK_MSG(DIM == 2 || DIM == 3, "quest::Triangle::area() only valid when dimension is 2 or 3");

      VectorType v(m_A, m_B);
      VectorType w(m_A, m_C);

      return (DIM==2)
                  ? 0.5 * std::fabs(v[0]*w[1] - v[1]*w[0])
                  : 0.5 * VectorType::cross_product( v,w).norm();
  }

  /*!
   * \brief Returns the barycentric coordinates of a triangle
   * \return The barycentric coordinates of the triangle inside a Point<T,3>
   */
  Point<T,3> barycenterCoords(const Point<T, DIM>& p) const{
	//adapted from Real Time Collision Detection by Christer Ericson

    VectorType u= VectorType(m_A, m_B);
    T nu, nv, ood;
    T x= std::abs(u[0]);
    T y= std::abs(u[1]);
    T z= std::abs(u[2]);


    if (x>=y && x>= z) {
        Point<T,2> pP= Point<T,2>::make_point(p[1], p[2]);
        Point<T,2> p1= Point<T,2>::make_point(m_A[1], m_A[2]);
        Point<T,2> p2= Point<T,2>::make_point(m_B[1], m_B[2]);
        Point<T,2> p3= Point<T,2>::make_point(m_C[1], m_C[2]);
        Triangle<T,2> t1(pP, p2, p3);
        Triangle<T,2> t2(pP, p3, p1);
        nu=t1.area2d();
        nv=t2.area2d();
        ood=1.0/u[0];

    }
    else if (y>=x && y>=z) {
        Point<T,2> pP= Point<T,2>::make_point(p[0], p[2]);
        Point<T,2> p1= Point<T,2>::make_point(m_A[0], m_A[2]);
        Point<T,2> p2= Point<T,2>::make_point(m_B[0], m_B[2]);
        Point<T,2> p3= Point<T,2>::make_point(m_C[0], m_C[2]);
        Triangle<T,2> t1(pP, p2, p3);
        Triangle<T,2> t2(pP, p3, p1);
        nu=t1.area2d();
        nv=t2.area2d();
        ood=-1.0/u[1];

    }
    else {
        Point<T,2> pP= Point<T,2>::make_point(p[0], p[1]);
        Point<T,2> p1= Point<T,2>::make_point(m_A[0], m_A[1]);
        Point<T,2> p2= Point<T,2>::make_point(m_B[0], m_B[1]);
        Point<T,2> p3= Point<T,2>::make_point(m_C[0], m_C[1]);
        Triangle<T,2> t1(pP, p2, p3);
        Triangle<T,2> t2(pP, p3, p1);
        nu=t1.area2d();
        nv=t2.area2d();
        ood=-1.0/u[2];
    }

    T bX= nu*ood;
    T bY = nv*ood;

    return (Point<T,3>::make_point(bX, bY, 1.0-bX-bY));
  }


  /*!
   *****************************************************************************
   * \brief Returns whether the triangle is degenerate
   * \return true iff the triangle is degenerate (0 area)
   * \see quest::Point
   *****************************************************************************
   */
  bool degenerate() const
  {
    return asctoolkit::utilities::isNearlyEqual(area(),  0.0, 1.0e-12);
  }


  /*!
   *****************************************************************************
   * \brief Returns whether Point P is in the triangle for some 3d Triangle
   * \return true iff P is in the triangle
   * \see quest::Point
   *****************************************************************************
   */

  bool checkInTriangle(const Point<double, DIM>& P) const{
    Point<T,3> bC= barycenterCoords(P);
    return ((bC[0]>=0.0) && (bC[1] >= 0.0) && (bC[2]>=0.0) &&
	    (bC[0]<=1.0) && (bC[1]<=1.0) && (bC[2]<=1.0) &&
	    (bC[0]+bC[1]+bC[2]<=1.0));
  }


  /*!
   * \brief Simple formatted print of a triangle instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
      os <<"{"
         << m_A <<" "
         << m_B <<" "
         << m_C <<"}";

      return os;
  }


private:

  PointType m_A;
  PointType m_B;
  PointType m_C;
};

//------------------------------------------------------------------------------
/// Free functions implementing Triangle's operators
//------------------------------------------------------------------------------

template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const Triangle<T,DIM> & tri)
{
    tri.print(os);
    return os;
}


} /* namespace quest */

#endif /* TRIANGLE_HPP_ */
