/*!
 * \file BezierCurve.hpp
 *
 * \brief A BezierCurve primitive for primal
 */

#ifndef PRIMAL_BEZIERCURVE_HPP_
#define PRIMAL_BEZIERCURVE_HPP_

#include "axom/slic.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/NumericArray.hpp"

#include "fmt/fmt.hpp"
#include <vector>
#include <ostream> // for std::ostream
#include <iostream>
#include <fstream>

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
  typedef std::vector< PointType > CoordsVec;

public:
  /*! Default constructor for an empty Bezier Curve*/
  BezierCurve() {}

  /*!
   * \brief Constructor for an empty Bezier Curve that reserves space for
   *  the given order of the curve
   *
   * \param [in] order the order of the resulting Bezier curve
   * \pre order is not negative
   *
   */
  BezierCurve(int ord)
  {
    SLIC_ASSERT(ord >= 0);
    m_controlpoints.reserve(ord+1);
    m_controlpoints.resize(ord+1);
  }
  
  /*!
   * \brief Constructor for a Bezier Curve from a list of Points 
   * \verbatim {x_0, x_1, x_2, x_3, y_0, y_1, y_2, y_3, z_0, z_1, z_2, z_3}  
   *
   * \param [in] pts an array with (n+1)*NDIMS entries, ordered by coordinate then by control point order
   * \param [in] ord number of control points minus 1 for which to reserve control point space
   * \pre order is not negative
   *
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
    for ( int i = 0 ; i <= ord ; i++)
    {
      for ( int j = 0 ; j < NDIMS ; j++)
      {
        tempar[j]=pts[j*(ord+1)+i];
      }
      this->addControlpoint(tempar);
    }
  }

  /*!
   * \brief Constructor for a Bezier Curve from an array of coordinates 
   *
   * \param [in] pts a vector with ord+1 points in it
   * \param [in] ord number of control points minus 1 for which to reserve control point space
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

    for (int i = 0 ; i <= ord ; i++)
    {
      this->addControlpoint(pts[i]);
    }
  }

  /*! Sets the order of the Bezier Curve*/
  void setOrder( int ord)
  { m_controlpoints.resize(ord+1); }

  /*! Returns the order of the Bezier Curve*/
  int getOrder() const
  { return m_controlpoints.size()-1; }

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

  /*! Retrieves the control point at index idx */
  PointType& operator[](int idx) { return m_controlpoints[idx]; }
  /*! Retrieves the control point at index idx */
  const PointType& operator[](int idx) const { return m_controlpoints[idx]; }
 
  std::vector< Point< T, NDIMS > > getControlPoints() const
  {
    return m_controlpoints; 
  }

  /*!
 * \brief Evaluates a Bezier Curve at particular parameter value between 0 and 1
 *
 * \param [in] t parameter value between 0 and 1 at which to evaluate
 * \return p the value of the Bezier Curve at t
 *
 */

  Point< T, NDIMS > eval_bezier(T t)
  {
    Point< T, NDIMS > ptval(0.0,NDIMS);
    int ord = m_controlpoints.size()-1;
    T dCarray[ord+1];
    
    for ( int i=0; i < NDIMS; i++)
    {
       for ( int j=0 ; j <= ord ; j++)
       {
         dCarray[j] = m_controlpoints[j][i];
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
  
  void split_bezier(T t, BezierCurve< T, NDIMS >& c1, BezierCurve< T, NDIMS>& c2)
  { 
    int ord = m_controlpoints.size()-1;
    c1.setOrder(ord);
    c2.setOrder(ord);
    T* dCarray = new T[NDIMS*2*(ord+1)];
    for ( int i=0; i < NDIMS; i++)
    {
       for ( int j=0 ; j <= ord ; j++)
       {
         dCarray[i*((ord+1))+j] = m_controlpoints[j][i];
       }
       dCarray[NDIMS*(ord+1)+(i)*((ord+1))]= dCarray[i*((ord+1))];
       for ( int j=1 ; j <= ord ; j++)
       {
         for ( int k=0 ; k <= ord-j ; k++)
         {
           dCarray[i*((ord+1))+k]=(1-t)*dCarray[i*((ord+1))+k]+t*dCarray[i*((ord+1))+k+1];
         }       
         dCarray[NDIMS*(ord+1)+(i)*((ord+1))+j]=dCarray[i*((ord+1))];
       }
    }
       c2=BezierCurve< T, NDIMS> (dCarray,ord);
       c1=BezierCurve< T, NDIMS> (dCarray+NDIMS*(ord+1),ord);
       delete [] dCarray;
       return;
  }

  /*!
   * \brief Simple formatted print of a Bezier Curve instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const 
  {
    const int sz = m_controlpoints.size()-1;

    os <<"{" << sz <<"-degree Bezier Curve:";
    for (int i=0 ; i< sz; ++i)
    {
      os << m_controlpoints[i] << ",";
    }
    os<<m_controlpoints[sz];
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
  CoordsVec m_controlpoints;
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
