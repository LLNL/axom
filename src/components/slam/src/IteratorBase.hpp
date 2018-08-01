/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file IteratorBase.hpp
 *
 * \brief Contains iterator base classes
 */

#ifndef SLAM_ITERBASE_H_
#define SLAM_ITERBASE_H_

#include "slam/Set.hpp"

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
 * It uses the Curiously recurring template pattern: The first template
 * parameter should be the actual iterator class. \n
 * The derived iterator class must implement `void advance(PositionType)`
 * to update the \a m_pos variable.
 * e.g.
 * \code{.cpp}
 * class IteratorDerived : public IteratorBase < IteratorDerived, DataType >
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
 * \tparam  DataType The return type of the iterator indirection
 */

template<class IterType, typename DataType>
class IteratorBase : public std::iterator<std::random_access_iterator_tag,
                                          DataType>
{
public:
  using PositionType = Set::PositionType;

protected:
  IteratorBase(int pos) : m_pos(pos) { }

private:
  /**
   *  Utility class to access protected function IterType::advance().
   *  Idea borrowed from: https://accu.org/index.php/journals/296
   */
  struct accessor : IterType
  {
    static void adv(IterType& derived, PositionType n)
    {
      void (IterType::* fn)(PositionType) = &accessor::advance;
      (derived.*fn)(n);
    }
  };

private:
  //call the derived iterator's advance() function
  void adv(IterType& derived, PositionType n) const {
    accessor::adv(derived, n);
  }

public:
  bool operator==(const IterType& other) const {
    return (m_pos == other.m_pos);
  }
  bool operator!=(const IterType& other) const { return !operator==(other); }
  bool operator<(const IterType& other) const { return m_pos < other.m_pos; }

  IterType& operator++() { adv(getIter(),1); return getIter(); }
  IterType operator++(int) {
    IterType ret = getIter();
    adv(getIter(), 1);
    return ret;
  }
  IterType& operator--() { adv(getIter(),-1); return getIter(); }
  IterType operator--(int) {
    IterType ret = getIter();
    adv(getIter(),-1);
    return ret;
  }

  IterType& operator+=(PositionType n) { adv(getIter(), n); return getIter(); }
  IterType& operator-=(PositionType n) { adv(getIter(), -n); return getIter(); }

  IterType operator+(PositionType n) const
  {
    IterType ret = getIter();
    adv(ret,n);
    return ret;
  }
  IterType operator-(PositionType n) const
  {
    IterType ret = getIter();
    adv(ret,-n);
    return ret;
  }

  friend PositionType operator-(const IterType& a, const IterType& b)
  {
    return (a.m_pos - b.m_pos);
  }

private:
  IterType& getIter() {
    return *static_cast<IterType*>(this);
  }
  const IterType& getIter() const {
    return *static_cast<const IterType*>(this);
  }

protected:
  int m_pos;
};

} // end namespace slam
} // end namespace axom

#endif //  SLAM_BITSET_H_
