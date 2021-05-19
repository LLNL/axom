// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IteratorBase.hpp
 *
 * \brief Contains iterator base classes
 */

#ifndef SLAM_ITERBASE_H_
#define SLAM_ITERBASE_H_

#include "axom/slam/Set.hpp"

namespace axom
{
namespace slam
{
/**
 * \class IteratorBase
 *
 * \brief Base class for a random access iterator over positions in a set
 *
 * This class is for keeping track of the position value in an iterator class.
 * It uses the Curiously Recurring Template Pattern (CRTP): The first template
 * parameter should be the actual iterator class. \n
 * The derived iterator class must implement `void advance(PositionType)`
 * to update the \a m_pos variable.
 * e.g.
 * \code{.cpp}
 * class IteratorDerived : public IteratorBase < IteratorDerived, PositionType >
 * {
 *   using IteratorBase<IteratorDerived, DataType>::m_pos;

 *   //...implementation of the accessing functions using m_pos...
 *
 *   //derive iterator class must implement this function, ideally protected.
 * protected:
 *   void advance( PositionType n) {...}
 * };
 * \endcode
 *
 * \tparam  IterType The iterator class to be the derived class
 * \tparam  PosType The type used to index into the iterator
 *          Must be a signed integral type
 */

template <typename IterType, typename PosType>
class IteratorBase
{
  static_assert(std::is_signed<PosType>::value, "PosType must be signed");
  static_assert(std::is_integral<PosType>::value, "PosType must be integral");

protected:
  IteratorBase() : m_pos(PosType()) { }
  explicit IteratorBase(PosType pos) : m_pos(pos) { }

private:
  /**
   *  Utility class to access protected function IterType::advance().
   *  Idea borrowed from: https://accu.org/index.php/journals/296
   */
  struct accessor : IterType
  {
    static void adv(IterType& derived, PosType n)
    {
      void (IterType::*fn)(PosType) = &accessor::advance;
      (derived.*fn)(n);
    }
  };

private:
  /// Call the derived iterator's advance() function
  void adv(IterType& derived, PosType n) const { accessor::adv(derived, n); }

public:
  /// \name Equality and relational operators
  /// \{

  /// Equality operator
  friend bool operator==(const IterType& lhs, const IterType& rhs)
  {
    return lhs.m_pos == rhs.m_pos;
  }
  /// Inequality operator
  friend bool operator!=(const IterType& lhs, const IterType& rhs)
  {
    return lhs.m_pos != rhs.m_pos;
  }
  /// Less than operator
  friend bool operator<(const IterType& lhs, const IterType& rhs)
  {
    return lhs.m_pos < rhs.m_pos;
  }
  /// Less than or equal operator
  friend bool operator<=(const IterType& lhs, const IterType& rhs)
  {
    return lhs.m_pos <= rhs.m_pos;
  }
  /// Greater than operator
  friend bool operator>(const IterType& lhs, const IterType& rhs)
  {
    return lhs.m_pos > rhs.m_pos;
  }
  /// Greater than or equal operator
  friend bool operator>=(const IterType& lhs, const IterType& rhs)
  {
    return lhs.m_pos >= rhs.m_pos;
  }
  /// \}

  /// \name Iterator advance and distance operators
  /// \{

  /// Pre-increment operator
  IterType& operator++()
  {
    adv(getIter(), 1);
    return getIter();
  }
  /// Post-increment operator
  IterType operator++(int)
  {
    IterType ret = getIter();
    adv(getIter(), 1);
    return ret;
  }
  /// Pre-decrement operator
  IterType& operator--()
  {
    adv(getIter(), -1);
    return getIter();
  }
  /// Post-decrement operator
  IterType operator--(int)
  {
    IterType ret = getIter();
    adv(getIter(), -1);
    return ret;
  }

  /// Addition-assignment operator
  IterType& operator+=(PosType n)
  {
    adv(getIter(), n);
    return getIter();
  }
  /// Subtraction-assignment operator
  IterType& operator-=(PosType n)
  {
    adv(getIter(), -n);
    return getIter();
  }

  /// Addition operator with iterator on left and position on right
  friend IterType operator+(const IterType& it, PosType n)
  {
    IterType ret(it);
    ret.adv(ret, n);
    return ret;
  }
  /// Addition operator with position on left and iterator on right
  friend IterType operator+(PosType n, const IterType& it)
  {
    return operator+(it, n);
  }

  /// Subtraction operator with iterator on left and position on right
  friend IterType operator-(const IterType& it, PosType n)
  {
    return operator+(it, -n);
  }
  /// Subtraction operator with position on left and iterator on right
  friend IterType operator-(PosType n, const IterType& it)
  {
    return operator+(it, -n);
  }

  /// Difference operator
  friend PosType operator-(const IterType& a, const IterType& b)
  {
    return (a.m_pos - b.m_pos);
  }
  /// \}

private:
  /// Accessor to derived class
  IterType& getIter() { return *static_cast<IterType*>(this); }
  /// Const accessor to derived class
  const IterType& getIter() const
  {
    return *static_cast<const IterType*>(this);
  }

protected:
  PosType m_pos;
};

}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_ITERATOR_BASE_H_
