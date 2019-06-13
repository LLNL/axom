// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef SLAM_POLICY_TRAITS_H_
#define SLAM_POLICY_TRAITS_H_

/**
 * \file PolicyTraits.hpp
 *
 * A collection of utility policy and traits classes for Slam.
 */

#include "axom/slam/Set.hpp"
#include "axom/slam/NullSet.hpp"
#include "axom/slam/policies/SizePolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

#include <utility>

namespace axom
{
namespace slam
{
namespace policies
{

/**
 * \brief Definition of a type trait to adapt a StridePolicy into a SizePolicy
 */
template<typename StridePolicyType, typename IntType, int VAL = 1>
struct StrideToSize
{
  using SizeType = CompileTimeSize<IntType, VAL>;
};

/**
 * \brief Specialization of StrideToSize trait for a RuntimeStride
 */
template<typename IntType>
struct StrideToSize < RuntimeStride<IntType>, IntType >
{
  using SizeType = RuntimeSize<IntType>;
};

/**
 * \brief Specialization of StrideToSize trait for a CompileTimeStride
 */
template<typename IntType, int VAL>
struct StrideToSize< CompileTimeStride<IntType, IntType(VAL)>, IntType, VAL >
{
  using SizeType = CompileTimeSize<IntType, IntType(VAL)>;
};

/**
 * \brief Specialization of StrideToSize trait for a StrideOne type
 */
template<typename IntType>
struct StrideToSize< StrideOne<IntType>, IntType >
{
  using SizeType = CompileTimeSize<IntType, StrideOne<IntType>::DEFAULT_VALUE >;
};


/**
 * \brief Type traits for null sets.
 *
 * The null pointer for most sets is nullptr
 */
template<typename SetType,
         typename P = typename SetType::PositionType,
         typename E = typename SetType::ElementType>
struct EmptySetTraits
{
  using EmptySetType = SetType;

  static EmptySetType* emptySet() { return nullptr; }

  static bool isEmpty(const EmptySetType* set)
  {
    return (set==emptySet() || set->empty() );
  }
};

/**
 * \brief Specialization of NullSetTraits for the base class Set
 *
 * The null pointer is of type NullSet
 */
template<typename P, typename E>
struct EmptySetTraits<slam::Set<P,E> >
{
  using EmptySetType = slam::Set<P,E>;

  static EmptySetType* emptySet()
  {
    static slam::NullSet<P,E> s_nullSet;
    return &s_nullSet;
  }

  static bool isEmpty(const EmptySetType* set)
  {
    return *set == *(emptySet()) || set->empty();
  }

};


} // end namespace policies

namespace traits
{
// Implementation of void_t from-C++17
template <class ... Ts>
using void_t = void;


///\name has_relation_ptr traits class
///@{

template<class T, class=void>
struct has_relation_ptr : std::false_type {};

template<class T>
struct has_relation_ptr<T, void_t<decltype(std::declval<T>().getRelation() )> >
  : std::true_type {};


///@}

///\name indices_use_indirection traits class for BivariateSetTypes
///@{

template<class T, class=void>
struct indices_use_indirection : std::true_type {};

template<class T>
struct indices_use_indirection<T, void_t<typename T::ProductSetType> >
  : std::false_type
{
  static_assert( std::is_base_of<typename T::BivariateSetType, T>::value, "" );
};


///@}

} // end namespace traits

} // end namespace slam
} // end namespace axom

#endif // SLAM_POLICY_TRAITS_H_
