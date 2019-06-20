/*!
 * \file BezierCurve.hpp
 *
 * \brief A BezierCurve primitive for primal
 */

#ifndef PRIMAL_BEZIERCURVE_HPP_
#define PRIMAL_BEZIERCURVE_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/NumericArray.hpp"

#include <vector>
#include <ostream> // for std::ostream

namespace axom
{
namespace primal
{ 

// Forward declare the templated classes and operator functions
template < typename T, int NDIMS >
class BezierCurve;

/*! \brief Overloaded output operator for Bezier Curves*/
template < typename T,int NDIMS >
std::ostream& operator<<(std::ostream & os, const BezierCurve< T,NDIMS > & bCurve);

/*!
 * \class BezierCurve 
 *
 * \brief Represents a Bezier Curve defined by an array of points
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 * \note The control points should be ordered from t=0 to t=1
 */

template < typename T,int NDIMS >
class BezierCurve
{
public:
  typedef Point< T,NDIMS >  PointType;
  typedef Vector< T,NDIMS > VectorType;
  typedef NumericArray< T,NDIMS > NumArrayType;

private:
  typedef std::vector< PointType > Coords;

public:
  /*! Default constructor for an empty Bezier Curve*/
  BezierCurve() {}

  /*!
   * \brief Constructor for an empty Bezier Curve that reserves space for
   *  the given order of the curve
   *
   * \param [in] order number of control points minus 1 for which to reserve space
   * \pre order is not negative
   *
   */
  BezierCurve(int order)
  {
    SLIC_ASSERT(order >= 0);
    m_controlpoints.reserve(order + 1);
  }

  /*! Return the order of the Bezier Curve*/
  int order() const { return m_controlpoints.size()-1; }

  /*! Appends a control point to the list of control points*/
  void addControlpoint(const PointType& pt)
  {
    m_controlpoints.push_back(pt);
  }

  /*! Clears the list of control points*/
  void clear()
  {
    m_controlpoints.clear();
  }

  /*! Retrieves the vertex at index idx */
  PointType& operator[](int idx) { return m_controlpoints[idx]; }
  /*! Retrieves the vertex at index idx */
  const PointType& operator[](int idx) const { return m_controlpoints[idx]; }

  /*!
   * \brief Simple formatted print of a Bezier Curve instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const 
  {
    const int sz = order();

    os <<"{" << sz <<"-order Bezier Curve:";
    for (int i=0 ; i< sz; ++i)
    {
      os << m_controlpoints[i] << ",";
    }
    if (sz >= 1)
    {
      os<<m_controlpoints[sz];
    }
    os<< "}";

    return os;
  }

  /*!
   * \brief Simple check for validity of a Bezier Curve
   *
   * Initial check is that the BezierCurve has two or more control points 
   * \return True, if the Bezier Curve is valid, False otherwise
   */
  bool isValid() const
  {
    return m_controlpoints.size() >= 2;
  }

private:
  Coords m_controlpoints;
}; // class BezierCurve

//------------------------------------------------------------------------------
/// Free functions implementing BezierCurve's operators
//------------------------------------------------------------------------------
template < typename T, int NDIMS >
std::ostream& operator<<(std::ostream & os, const BezierCurve< T,NDIMS > & bCurve)
{
  bCurve.print(os);
  return os;
} // std::ostream& operator <<...

} // namespace primal
} // namespace axom

#endif // PRIMAL_BEZIERCURVE_HPP_
