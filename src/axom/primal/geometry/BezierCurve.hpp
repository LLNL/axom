// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file BezierCurve.hpp
 *
 * \brief A BezierCurve primitive
 */

#ifndef PRIMAL_BEZIERCURVE_HPP_
#define PRIMAL_BEZIERCURVE_HPP_

#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "fmt/fmt.hpp"
#include <vector>
#include <ostream>

namespace axom
{
namespace primal
{

// Forward declare the templated classes and operator functions
template < typename T, int NDIMS >
class BezierCurve;

/*! \brief Overloaded output operator for Bezier Curves*/
template < typename T,int NDIMS >
std::ostream& operator<<(std::ostream & os,
                         const BezierCurve< T,NDIMS > & bCurve);

// UedaAreaMats is where all of the area matrices from Ueda '99 are stored
extern std::map<int, std::vector<double> > UedaAreaMats;

/*!
 * \class BezierCurve
 *
 * \brief Represents a Bezier curve defined by an array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * \note The order of a Bezier curve with N+1 control points is N
 * \note The control points should be ordered from t=0 to t=1
 */

template < typename T,int NDIMS >
class BezierCurve
{
public:
  using PointType = Point< T,NDIMS >;
  using VectorType = Vector< T,NDIMS >;
  using NumArrayType = NumericArray< T,NDIMS >;
  using SegmentType = Segment< T, NDIMS >;
  using CoordsVec = std::vector< PointType >;

public:
  /*! Default constructor for an empty Bezier Curve*/
  BezierCurve() = default;

  /*!
   * \brief Constructor for an empty Bezier Curve that reserves space for
   *  the given order of the curve
   *
   * \param [in] order the order of the resulting Bezier curve
   * \pre order is not negative
   */
  BezierCurve(int ord)
  {
    SLIC_ASSERT(ord >= 0);
    m_controlpoints.reserve(ord+1);
    m_controlpoints.resize(ord+1);
  }

  /*!
   * \brief Constructor for a Bezier Curve from a list of Points
   * \verbatim {x_0, x_1, x_2, x_3,
   *            y_0, y_1, y_2, y_3,
   *            z_0, z_1, z_2, z_3}
   *
   * \param [in] pts an array with (n+1)*NDIMS entries, ordered by coordinate
   * then by control point order
   * \param [in] ord number of control points minus 1 for which to reserve
   * control point space
   * \pre order is not negative
   */
  BezierCurve(T* pts, int ord)
  {
    if ( ord <= 0 )
    {
      clear();
    }
    // sanity check
    SLIC_ASSERT(pts != nullptr);

    m_controlpoints.reserve(ord+1);

    T tempar[NDIMS];
    for ( int p = 0 ; p <= ord ; p++)
    {
      for ( int j = 0 ; j < NDIMS ; j++)
      {
        tempar[j]=pts[j*(ord+1)+p];
      }
      this->addControlPoint(tempar);
    }
  }

  /*!
   * \brief Constructor for a Bezier Curve from an array of coordinates
   *
   * \param [in] pts a vector with ord+1 points in it
   * \param [in] ord number of control points minus 1 for which to reserve
   * control point space
   * \pre order is not negative
   *
   */

  BezierCurve(PointType* pts, int ord)
  {
    if (ord <= 0)
    {
      clear();
    }
    // sanity check
    SLIC_ASSERT(pts != nullptr);

    m_controlpoints.reserve(ord+1);

    for (int p = 0 ; p <= ord ; ++p)
    {
      this->addControlPoint(pts[p]);
    }
  }

  /*! Sets the order of the Bezier Curve*/
  void setOrder( int ord)
  { m_controlpoints.resize(ord+1); }

  /*! Returns the order of the Bezier Curve*/
  int getOrder() const
  { return m_controlpoints.size()-1; }

  /*! Appends a control point to the list of control points*/
  void addControlPoint(const PointType& pt)
  {
    m_controlpoints.push_back(pt);
  }

  /*! Clears the list of control points*/
  void clear()
  {
    m_controlpoints.clear();
  }

  /*! Retrieves the control point at index \a idx */
  PointType& operator[](int idx) { return m_controlpoints[idx]; }
  /*! Retrieves the control point at index \a idx */
  const PointType& operator[](int idx) const { return m_controlpoints[idx]; }

  /* Checks equality of two Bezier Curve */
  friend inline bool operator==(const BezierCurve< T, NDIMS>& lhs,
                                const BezierCurve< T, NDIMS>& rhs)
  {
    return lhs.m_controlpoints == rhs.m_controlpoints;
  }

  friend inline bool operator!=(const BezierCurve< T, NDIMS>& lhs,
                                const BezierCurve< T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  CoordsVec getControlPoints() const
  {
    return m_controlpoints;
  }

  /*!
   * \brief Evaluates a Bezier curve at a particular parameter value \a t
   *
   * \param [in] t parameter value at which to evaluate
   * \return p the value of the Bezier curve at t
   *
   * \note We typically evaluate the curve at \a t between 0 and 1
   */

  PointType evaluate(T t) const
  {
    PointType ptval;

    const int ord = getOrder();
    std::vector<T> dCarray(ord+1);

    // Run de Casteljau algorithm on each dimension
    for ( int i=0 ; i < NDIMS ; ++i)
    {
      for ( int p=0 ; p <= ord ; ++p)
      {
        dCarray[p] = m_controlpoints[p][i];
      }

      for ( int p=1 ; p <= ord ; ++p)
      {
        const int end = ord-p;
        for ( int k=0 ; k <= end ; ++k)
        {
          dCarray[k]=(1-t)*dCarray[k] + t*dCarray[k+1];
        }
      }
      ptval[i]=dCarray[0];
    }

    return ptval;
  }

  /*!
   * \brief Splits a Bezier curve into two Bezier curves at particular parameter
   * value between 0 and 1
   *
   * \param [in] t parameter value between 0 and 1 at which to evaluate
   * \param [out] c1, c2 Bezier curves that split the original
   */
  void split(T t, BezierCurve& c1, BezierCurve& c2) const
  {
    int ord = getOrder();
    SLIC_ASSERT(ord >= 0);

    // Note: the second curve's control points are computed inline
    //       as we find the first curve's control points
    c2 = *this;

    c1.setOrder(ord);
    c1[0] = c2[0];

    // Run de Casteljau algorithm
    // After each iteration, save the first control point into c1
    for ( int p=1 ; p <= ord ; ++p)
    {
      const int end = ord-p;
      for(int k=0 ; k<= end ; ++k)
      {
        PointType& pt1 = c2[k];
        const PointType& pt2 = c2[k+1];
        for(int i=0 ; i< NDIMS ; ++i)
        {
          pt1[i] = (1-t)*pt1[i] + t*pt2[i];
        }
      }
      c1[p] = c2[0];
    }

    return;
  }


 /* 
   * \brief Calculates the sector area (area between curve and origin) of a Bezier Curve
   *
   * \param [in] tol a tolerance parameter controlling definition of
   * near-linearity
   * \param [out] boolean TRUE if c1 is near-linear
   */
  T sectorArea() const
  {
    T A = 0;
        int ord = getOrder();
        if (UedaAreaMats.find(ord)==UedaAreaMats.end())
        {
          std::vector<double>  newUedaAreaMat((ord+1)*(ord+1));
          int twonchoosen=axom::utilities::binomial_coefficient(2*ord,ord);
          for (int i=0; i<=ord; ++i)
          { 
            for (int j=0; j<=ord; ++j)
            {
              if ((i==0 && j==0) || (i==(ord)&& j==(ord)))
              {newUedaAreaMat[i*(ord+1)+j]=0.0;}
              else
              {
              newUedaAreaMat[i*(ord+1)+j]=((1.0*j-i)/2)*(2.0*(ord)/(1.0*twonchoosen))*
                (1.0*axom::utilities::binomial_coefficient(i+j,i)/(1.0*i+j))*
                (1.0*axom::utilities::binomial_coefficient(2*(ord)-i-j,(ord)-j)/(1.0*(ord)-j+(ord)-i));
              }
            } 
          }
          UedaAreaMats.insert(std::pair<int,std::vector<double>>(ord,newUedaAreaMat));
        }
        const std::vector<double> &whicharea = (UedaAreaMats.find(ord)->second);
        for (int i=0; i<=ord; ++i)
        { 
          for (int j=0; j<=ord; ++j)
          {
            A+=static_cast<T>(whicharea[i*(ord+1)+j])*m_controlpoints[i][1]*m_controlpoints[j][0];
          } 
        }
      return A;
  } 
  
  /*!
   * \brief Predicate to check if the Bezier curve is approximately linear
   *
   * \param [in] tol a tolerance parameter controlling definition of
   * near-linearity
   * \param [out] boolean TRUE if c1 is near-linear
   */

  bool is_linear(double tol) const
  {
    int ord1 = m_controlpoints.size()-1;
    SegmentType linear1(m_controlpoints[0],m_controlpoints[ord1]);
    double d1=0.0;
    for (int p=1 ; p<ord1 ; ++p) // line defined by endpoints; only add interior
    {
      d1=d1+squared_distance(m_controlpoints[p],linear1);
    }
    return (d1<tol); // Note: tolerance is squared
  }

  /*!
   * \brief Simple formatted print of a Bezier Curve instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */

  std::ostream& print(std::ostream& os) const
  {
    const int ord = getOrder();

    os <<"{" << ord <<"-degree Bezier Curve:";
    for (int p=0 ; p< ord ; ++p)
    {
      os << m_controlpoints[p] << ",";
    }
    os<<m_controlpoints[ord];
    os<< "}";

    return os;
  }

private:
  CoordsVec m_controlpoints;
}; // class BezierCurve

//------------------------------------------------------------------------------
/// Free functions implementing BezierCurve's operators
//------------------------------------------------------------------------------
template < typename T, int NDIMS >
std::ostream& operator<<(std::ostream & os,
                         const BezierCurve< T,NDIMS > & bCurve)
{
  bCurve.print(os);
  return os;
}

} // namespace primal
} // namespace axom

#endif // PRIMAL_BEZIERCURVE_HPP_
