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

namespace axom
{
namespace slam
{

/**
 * \class IteratorBase
 *
 * \brief An iterator base class that keeps track of the position value.
 *
 * This class is for keeping track of the position value in an iterator class.
 * It uses the Curiously recurring template pattern: The first template
 * parameter should be the actual iterator class. \n
 * e.g.
 * \code{.cpp}
 * class IteratorDerived : public IteratorBase < IteratorDerived, DataType >
 * {
 *   using IteratorBase<IteratorDerived, DataType>::m_pos;

 *   //...implementation of the accessing functions using m_pos...
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
  using PositionType = MeshIndexType;

protected:
  IteratorBase(int pos) : m_pos(pos) { }

public:
  bool operator==(const IterType& other) const {
    return (m_pos == other.m_pos);
  }
  bool operator!=(const IterType& other) const { return !operator==(other); }
  bool operator<(const IterType& other) const { return m_pos < other.m_pos; }

  IterType& operator++() { getRetIter().advance(1); return getRetIter(); }
  IterType operator++(int) {
    IterType ret = getRetIter();
    getRetIter().advance(1); return ret;
  }
  IterType& operator--() { getRetIter().advance(-1); return getRetIter(); }
  IterType operator--(int) {
    IterType ret = getRetIter(); getRetIter().advance(
      -1); return ret;
  }

  IterType& operator+=(PositionType n) {
    getRetIter().advance(n);
    return getRetIter();
  }
  IterType& operator-=(PositionType n) {
    getRetIter().advance(-n);
    return getRetIter();
  }

  IterType operator+(PositionType n) const
  {
    IterType ret = getRetIter(); ret.advance(n);
    return ret;
  }
  IterType operator-(PositionType n) const
  {
    IterType ret = getRetIter(); ret.advance(-n);
    return ret;
  }

  friend PositionType operator-(const IterType& a, const IterType& b)
  {
    return (a.m_pos - b.m_pos);
  }

private:
  IterType& getRetIter() { return *((IterType*)this); }
  const IterType& getRetIter() const { return *((IterType*)this); }

protected:
  /** Derived iterator class can override this function for specialized
   * implementation. (But the function has to be declared public to compile)
   */
  virtual void advance(PositionType n) { m_pos += n; }

protected:
  int m_pos;
};

} // end namespace slam
} // end namespace axom

#endif //  SLAM_BITSET_H_
