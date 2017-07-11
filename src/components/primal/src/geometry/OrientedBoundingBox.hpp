/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef ORIENTEDBOUNDINGBOX_HPP_
#define ORIENTEDBOUNDINGBOX_HPP_

#include <vector>

#include "axom/config.hpp"   // defines AXOM_USE_CXX11

#include "primal/NumericArray.hpp" // for numeric arrays
#include "axom_utils/Utilities.hpp" // for nearly equal
#include "axom_utils/Matrix.hpp" // for Matrix
#include "axom_utils/eigen_solve.hpp" // for eigen_solve
#include "primal/Point.hpp"
#include "primal/Vector.hpp"
#include "primal/BoundingBox.hpp"

#include "slic/slic.hpp"

namespace axom {
namespace primal {

// Forward declare the templated classes and operator functions
template < typename T,int NDIMS >
class OrientedBoundingBox;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 *******************************************************************************
 * \brief Equality comparison operator for oriented bounding boxes.
 * Two oriented bounding boxes are equal when they have the same bounds and axes.
 *******************************************************************************
 */
template < typename T,int NDIMS >
bool operator==( const OrientedBoundingBox< T, NDIMS > &lhs,
                 const OrientedBoundingBox< T, NDIMS > &rhs );

/*!
 *******************************************************************************
 * \brief Inequality comparison operator for oriented bounding boxes.
 * Two bounding boxes are unequal when they have different bounds
 *******************************************************************************
 */
template < typename T,int NDIMS >
bool operator!=( const OrientedBoundingBox< T, NDIMS > & lhs,
                 const OrientedBoundingBox< T, NDIMS >& rhs   );

/*!
 *******************************************************************************
 * \brief Overloaded output operator for oriented bounding boxes
 *******************************************************************************
 */
template < typename T,int NDIMS >
std::ostream& operator<<( std::ostream & os,
                          const OrientedBoundingBox< T,NDIMS >& box );
///@}

/*!
 *******************************************************************************
 * \class
 *
 * \brief OrientedBoundingBox represents an oriented bounding box defined by
 * its min and max coordinates and it's axes.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *******************************************************************************
 */

template < typename T, int NDIMS >
class OrientedBoundingBox
{
public:
  typedef Point< T, NDIMS > PointType;
  typedef Vector< T, NDIMS > VectorType;
  typedef OrientedBoundingBox< T, NDIMS > OrientedBoxType;
  typedef BoundingBox< T, NDIMS > BoxType;

public:

  /*!
   *****************************************************************************
   * \brief Constructor. Creates an invalid oriented bounding box by setting
   * the first coordinate of the first axis to the smallest number of type T.
   ****************************************************************************
   */
  OrientedBoundingBox();

  /*!
   *****************************************************************************
   * \brief Constructor. Creates an oriented bounding box containing a single
   * point. Extents are all zero, and axes are set to identity.
   *****************************************************************************
   */
  OrientedBoundingBox( const PointType& pt );

  /*!
   *****************************************************************************
   * \brief Constructor. Creates an oriented bounding box with given center,
   * axes, and extents. Normalizes axes and sets any negative extent to 0.
   * \param [in] c center of OBB
   * \param [in] u axes of OBB
   * \param [in] e extents of OBB
   * \note Axes are made orthonormal, but if they don't span, some will be
   * set to 0.
   *****************************************************************************
   */
  OrientedBoundingBox( const PointType &c, const VectorType u[NDIMS],
                       const VectorType &e );

  /*!
   *****************************************************************************
   * \brief Copy Constructor.
   * \param [in] other The oriented bounding box to copy
   *****************************************************************************
   */
  OrientedBoundingBox( const OrientedBoundingBox& other );

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  ~OrientedBoundingBox() { }

  /*!
   *****************************************************************************
   * \brief Resets the bounds to those of the default constructor
   * \note This invalidates the bounding box (i.e. isValid() will be false)
   *****************************************************************************
   */
  void clear();

  /*!
   *****************************************************************************
   * \brief Returns the centroid (midpoint) of the bounding box.
   * \return Point at the bounding box centroid.
   *****************************************************************************
   */
  PointType centroid() const { return m_c; }

  /*!
   *****************************************************************************
   * \brief Return the axes
   * \param [in] u. Array of vectors the result will be passed through.
   *****************************************************************************
   */
  void axes(VectorType u[NDIMS]) const;

  /*!
   *****************************************************************************
   * \brief Returns the extents of the oriented bounding box.
   * \return extents in a vector
   *****************************************************************************
   */
  VectorType extents() const { return m_e; }

  /*!
   *****************************************************************************
   * \brief Returns the vertices of the oriented bounding box.
   * \return vertices in a std::vector<VectorType>
   *****************************************************************************
   */
  std::vector< PointType > vertices() const;

  
  /*!
   *****************************************************************************
   * \brief Returns the dimension of the ambient space for this bounding box.
   * \return d the dimension of this bounding box instance.
   * \post d >= 1.
   *****************************************************************************
   */
  int dimension() const { return NDIMS; };

  /*!
   *****************************************************************************
   * \brief Expands the box so it contains the passed in point.
   * \param pt [in] point to add in
   * \post this->contains(pt) == true
   *****************************************************************************
   */
  void addPoint(PointType pt);

  /*!
   *****************************************************************************
   * \brief Expands the box so it contains the passed in box.
   * \param obb [in] OBB to add in
   * \note If obb is invalid, makes this invalid too.
   * \post this->contains(obb) == true
   *****************************************************************************
   */
  void addBox(OrientedBoxType obb);

  /*!
   *****************************************************************************
   * \brief Expands the extents by the given amount.
   * \param [in] expansionAmount an absolute amount to expand
   * \note This function checks to ensure the bounding box is valid afterwards.
   * \return A reference to the bounding box after it has been expanded
   *****************************************************************************
   */
  OrientedBoundingBox& expand(T expansionAmount);

  /*!
   *****************************************************************************
   * \brief Scales the bounding box about its center by a given amount.
   * \param [in] scaleFactor the multiplicative factor by which to scale
   * \note Checks to ensure that the bounding box is valid after inflation.
   * \note If scaleFactor is less than 1, the bounding box will shrink.
   * \note If scaleFactor is 0, the bounding box will shrink to its midpoint
   *  towards the center, and we fix the bounds after shrinking
   * \note if scaleFactor is less than 0, we use -scaleFactor.
   *  \return A reference to the bounding box after it has been scaled
   *****************************************************************************
   */
  OrientedBoundingBox& scale(double scaleFactor);

  /*!
   *****************************************************************************
   * \brief Shifts the bounding box by a fixed displacement.
   * \param [in] displacement the amount with which to move the bounding box
   * \note This function checks to ensure the bounding box is valid afterwards.
   * \return A reference to the bounding box after it has been shifted
   *****************************************************************************
   */
  OrientedBoundingBox& shift(const VectorType& displacement);

  /*!
   *****************************************************************************
   * \brief Overloaded assignment operator.
   * \param [in] rhs oriented bounding box instance on the right-hand side
   * \return
   *****************************************************************************
   */
  OrientedBoundingBox& operator=(const OrientedBoundingBox& rhs );

  /*!
   *****************************************************************************
   * \brief Checks whether the box contains the point
   * \param [in] otherPt the point that we are checking
   * \return status true if point inside the box, else false.
   *
   * \note This function assumes all intervals are closed
   * (i.e. contain their boundaries).  We may need to deal with open
   * and half open boundaries in the future.
   *****************************************************************************
   */
  template < typename OtherType >
  bool contains( const Point< OtherType, NDIMS >& otherPt, double EPS=1E-4)
    const;

  /*!
   *****************************************************************************
   * \brief Checks whether the box fully contains another bounding box
   * \param [in] otherOBB the bounding box that we are checking
   * \return status true if bb is inside the box, else false.
   * \note We are allowing the other bounding box to have a different coordinate
   *  type. This should work as long as the two Ts are comparable with
   *  operator<().
   *****************************************************************************
   */
  template < typename OtherType >
  bool contains( const OrientedBoundingBox< OtherType, NDIMS >& otherOBB ) const;

  /*!
   *****************************************************************************
   * \brief Returns the same point but in the local coordinates of this.
   * \param [in] pt the pt we want to find in local coordinates
   * \return Vector which is just pt in local coordinates, i.e. is the
   *  projection of the difference between pt and m_c onto the axes m_u.
   *****************************************************************************
   */
  VectorType toLocal( const PointType &pt ) const;

  /*!
   *****************************************************************************
   * \brief Returns the furthest point in this OBB from the passed in pt.
   * \param [in] pt the pt we want to find the furthest point from
   * \note If there are multiple furthest point the are no guarantees
   * on which one it returns..
   * \return furthest point in this OBB from pt
   *****************************************************************************
   */
  PointType furthestPoint( const PointType &pt ) const;

  /*!
   *****************************************************************************
   * \brief Checks that we have a valid bounding box.
   *****************************************************************************
   */
  bool isValid() const;

  /*!
   *****************************************************************************
   * \brief Subdivides this oriented bounding box instance into two sub-boxes by
   *  splitting along the longest dimension.
   * \param [in,out] right the right sub-box.
   * \param [in,out] left  the left sub-box.
   *****************************************************************************
   */
  void bisect( OrientedBoxType& right, OrientedBoxType& left ) const;

  /*!
   *****************************************************************************
   * \brief Simple formatted print of a bounding box instance.
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   *****************************************************************************
   */
  std::ostream& print( std::ostream& os ) const;

  /// @}

 private:

  /// \name Internal Helper Routines
  /// @{

  /*!
   *****************************************************************************
   * \brief Recursively enumerates the vertices of this, storing them in l.
   *****************************************************************************
   */ 
  static void vertex_enum(std::vector< Point< T, NDIMS > > &l, int i,
    const Vector< T, NDIMS > u[NDIMS], const Vector< T, NDIMS > &e,
    Vector< T, NDIMS > curr)
  {

    if (i == NDIMS) {  // base case
      l.push_back(Point< T, NDIMS >(curr.array()));
    } else {
      curr += e[i]*u[i];
      vertex_enum(l, i + 1, u, e, curr);
      curr -= (static_cast< T >(2.)*e[i])*u[i];
      vertex_enum(l, i + 1, u, e, curr);
    }
  }

  /// @}

  /*!
   *****************************************************************************
   * \brief Ensures the axes are an orthonormal basis and all extents are
   * nonnegative. Applies Gram-Schmidt to them. If any extent is negative it
   * sets it to 0.
   *****************************************************************************
   */
  void checkAndFix();


 private:
  PointType m_c;  // centroid
  VectorType m_u[NDIMS];  // axes (all unit length)
  VectorType m_e;  // extents (each coord >= 0 or it's invalid)
};

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  OrientedBoundingBox implementation
//------------------------------------------------------------------------------

namespace axom {
namespace primal {


//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >::OrientedBoundingBox()
{
  (this->m_u[0])[0] = ValueRange< T >::lowest();
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >::OrientedBoundingBox(
  const Point< T, NDIMS > &pt )
{
  this->m_c = Point< T, NDIMS >(pt);

  for (int i = 0; i < NDIMS; i++) {
    // initialize ith axis to ith standard basis vector
    this->m_u[i] = Vector< T, NDIMS >();
    this->m_u[i][i] = static_cast< T >(1.);

    // initialize ith extent to 0
    this->m_e[i] = T();
  }
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >::OrientedBoundingBox( const Point< T,NDIMS > &c,
  const Vector< T,NDIMS > u[NDIMS], const Vector< T,NDIMS > &e )
{
  this->m_c = Point< T, NDIMS >(c);
  for (int i = 0; i < NDIMS; i++) this->m_u[i] = Vector< T, NDIMS >(u[i]);
  this->m_e = Vector< T, NDIMS >(e);
  this->checkAndFix();
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >::OrientedBoundingBox(
  const OrientedBoundingBox& other )
{
  this->m_c = other.centroid();
  this->m_e = other.extents();
  other.axes(this->m_u);
  this->checkAndFix();
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
void OrientedBoundingBox< T, NDIMS >::clear()
{
  (this->m_u[0])[0] = ValueRange< T >::lowest();
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
void OrientedBoundingBox< T, NDIMS >::axes( Vector< T, NDIMS > u[NDIMS] ) const
{
  for (int i = 0; i < NDIMS; i++) u[i] = this->m_u[i];
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
void OrientedBoundingBox< T, NDIMS >::addPoint( Point< T, NDIMS > pt )
{
  Vector< T, NDIMS > pt_l = this->toLocal(pt);
  for (int i = 0; i < NDIMS; i++) {
    T proj = pt_l[i];
    if (proj < T()) proj = -proj;

    if (proj > this->m_e[i]) this->m_e[i] = proj;
  }
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
void OrientedBoundingBox< T, NDIMS >::addBox( OrientedBoundingBox< T, NDIMS > obb )
{
  SLIC_CHECK_MSG(obb.isValid(), "Passed in OBB is invalid.");
  if (!obb.isValid()) {
    // make this invalid
    this->clear();
    return;
  }

  std::vector< Point< T, NDIMS > > res = obb.vertices();
  int size = res.size();

  for (int i = 0; i < size; i++)
    this->addPoint(res[i]);

}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
std::vector< Point< T, NDIMS > > OrientedBoundingBox< T, NDIMS >::vertices()
  const
{
  std::vector< Point< T, NDIMS > > res;
  Vector< T, NDIMS > curr(this->m_c);

  OrientedBoundingBox< T, NDIMS >::vertex_enum(res, 0, this->m_u, this->m_e,
    curr);
  return res;
}


//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >& OrientedBoundingBox< T, NDIMS >::expand(
  T expansionAmount)
{
  for (int i = 0; i < NDIMS; i++) this->m_e[i] += expansionAmount;
  this->checkAndFix();

  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >& OrientedBoundingBox< T, NDIMS >::scale(
  double scaleFactor)
{
  if (scaleFactor < T()) scaleFactor = -scaleFactor;
  for (int i = 0; i < NDIMS; i++) this->m_e[i] *= scaleFactor;
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >& OrientedBoundingBox< T, NDIMS >::shift(
  const Vector< T, NDIMS >& displacement)
{
  NumericArray< T, NDIMS > &coords = (this->m_c).array();
  NumericArray< T, NDIMS > c;
  for (int i = 0; i < NDIMS; i++) c[i] = coords[i] + displacement[i];
  this->m_c = Point< T, NDIMS >(c);
  checkAndFix();
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS >& OrientedBoundingBox< T, NDIMS >::operator=(
  const OrientedBoundingBox& rhs )
{
  if (this != &rhs) {
    this->m_e = rhs.extents();
    this->m_c = rhs.centroid();
    rhs.axes(this->m_u);
    this->checkAndFix();
  }
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
template < typename OtherType >
bool OrientedBoundingBox< T, NDIMS >::contains(
  const Point< OtherType, NDIMS >& otherPt, double EPS) const
{
  if (!(this->isValid())) return false;
  Vector< T, NDIMS > pt_l = this->toLocal(otherPt);
  T proj;
  T margin = static_cast< T >(EPS);
  for (int i = 0; i < NDIMS; i++) {
    proj = pt_l[i];

    // make proj nonnegative
    if (proj < T()) proj = -proj;

    if (proj > (m_e[i] + margin)) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
template < typename OtherType >
bool OrientedBoundingBox< T, NDIMS >::contains(
  const OrientedBoundingBox< OtherType, NDIMS >& otherOBB) const
{
  std::vector< Point< T, NDIMS > > l = otherOBB.vertices();

  int size = l.size();

  for (int i = 0; i < size; i++) {
    if (!this->contains(l[i])) return false;
  }

  return true;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Vector< T, NDIMS > OrientedBoundingBox< T, NDIMS >::toLocal(
  const Point< T, NDIMS > &pt ) const
{
  Vector< T, NDIMS > d(pt);
  for (int i = 0; i < NDIMS; i++) {
    d[i] -= (this->m_c[i]);
  }

  Vector< T, NDIMS > res;
  for (int i = 0; i < NDIMS; i++) {
    res[i] = d.dot(this->m_u[i]);
  }
  return res;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Point< T, NDIMS > OrientedBoundingBox< T, NDIMS >::furthestPoint(
  const PointType &pt ) const
{
  // TODO: toLocal()
  Vector< T, NDIMS > d = Vector< T, NDIMS >(pt) - Vector<T, NDIMS >(this->m_c);
  Vector< T, NDIMS > res(this->m_c);
  for (int i = 0; i < NDIMS; i++) {
    // since the local coordinates are individually constrained, we can simply
    // choose maximize the objective in each direction independently, meaning
    // meaning there is a simple analytical solution
    T dot = d.dot(this->m_u[i]);

    if (dot > T()) {
      res -= (this->m_e[i])*(this->m_u[i]);
    } else {
      res += (this->m_e[i])*(this->m_u[i]);
    }
  }
  return Point< T, NDIMS >(res.array());
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
bool OrientedBoundingBox< T, NDIMS >::isValid() const
{
  return !((this->m_u[0])[0] == ValueRange< T >::lowest());
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
void OrientedBoundingBox< T, NDIMS >::bisect(
  OrientedBoxType& right, OrientedBoxType& left) const
{
  T max = this->m_e[0];
  int bi = 0;
  for (int i = 1; i < NDIMS; i++) {
    if (max < this->m_e[i]) {
      max = this->m_e[i];
      bi = i;
    }
  }


  Vector< T, NDIMS > c1 = Vector< T, NDIMS >(this->m_c);
  Vector< T, NDIMS > c2 = Vector< T, NDIMS >(this->m_c);

  c1 += (static_cast< T >(0.5)*this->m_e[bi])*m_u[bi];
  c2 -= (static_cast< T >(0.5)*this->m_e[bi])*m_u[bi];

  Vector< T, NDIMS > e(this->m_e);
  e[bi] *= (static_cast< T >(0.5));

  right = OrientedBoundingBox< T, NDIMS >(Point< T, NDIMS >(c1.array()),
    this->m_u, e);
  left = OrientedBoundingBox< T, NDIMS >(Point< T, NDIMS >(c2.array()),
    this->m_u, e);
}


//------------------------------------------------------------------------------
template < typename T, int NDIMS >
std::ostream& OrientedBoundingBox< T, NDIMS >::print( std::ostream& os ) const
{
  // TODO: see if I can improve this
  os << "{centroid:" << this->m_c << "; extents:" << this->m_e;
  for (int i = 0; i < NDIMS; i++) {
    os << " axis " << i << ":" << this->m_u[i];
  }
  os << "}";
  return os;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
void OrientedBoundingBox< T, NDIMS >::checkAndFix()
{
  // if this is invalid, don't do anything!
  if (!this->isValid()) return;

  // TODO: talk with someone about whether this EPS thing makes sense.
  // Got the idea from Vector.hpp
  static const double EPS = 1.0e-8;

  // do Gram-Schmidt
  for (int i = 0; i < NDIMS; i++) {
    Vector< T, NDIMS > temp(this->m_u[i]);

    for (int j = 0; j < i; j++)
      temp -= (this->m_u[i].dot(this->m_u[j]))*(this->m_u[j]);

    double norm = temp.norm();

    if (norm < EPS) {
      // this one lies in the span of the others
      // so set it to 0
      this->m_u[i] = Vector< T, NDIMS >();
      continue;
    } else if (norm - 1. > EPS || norm - 1. < -EPS) {
      for (int k = 0; k < NDIMS; k++) {
        this->m_u[i][k] = static_cast< T >(
          static_cast< double >(temp[k])/norm);
      }
    }
  }

  // check and fix extents
  for (int i = 0; i < NDIMS; i++)
    if (this->m_e[i] < T()) this->m_e[i] = T();
}


//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template < typename T, int NDIMS >
bool operator==( const OrientedBoundingBox< T, NDIMS > &lhs,
                 const OrientedBoundingBox< T, NDIMS > &rhs )
{
  return lhs.contains(rhs) && rhs.contains(lhs);
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
bool operator!=( const OrientedBoundingBox< T, NDIMS > & lhs,
                 const OrientedBoundingBox< T, NDIMS >& rhs   )
{
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
std::ostream& operator<<( std::ostream & os,
                          const OrientedBoundingBox< T,NDIMS >& box )
{
  return box.print(os);
}



} /* namespace primal */
} /* namespace axom */

#endif /* ORIENTEDBOUNDINGBOX_HPP_ */






