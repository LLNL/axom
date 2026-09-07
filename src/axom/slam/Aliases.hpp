// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Aliases.hpp
 * \brief Common static-relation configurations over external buffers.
 *
 * VariableRelation, ConstantRelation, and RuntimeConstantRelation bind axom::Array
 * objects by pointer. Their View counterparts store axom::ArrayView values.
 * Neither form owns its buffers or its from/to sets. Array objects must outlive
 * Array-bound relations and allocations must outlive view-bound relations.
 * Referenced sets must outlive either form.
 *
 * Relation entries are ToSet::PositionType, rather than ToSet::ElementType.
 * FlatPosType indexes the relation's flat storage and its begin offsets.
 * Its default can represent both sets' positions, so choose a wider signed type
 * when the total number of relation entries needs it.
 *
 * View names select a buffer representation. Referenced sets and buffers 
 * must be accessible to the code using the relation.
 *
 * Set aliases live in RangeSet.hpp and IndirectionSet.hpp.
 * We do not provide map aliases here since Map and BivariateMap already default
 * to Array storage and stride one, and more complex configurations using
 * STLVectorIndirection, CArrayIndirection, or MappedVariableCardinality
 * are just as easy to represent directly.
 */

#pragma once

#include "axom/slam/StaticRelation.hpp"

#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/StridePolicies.hpp"

namespace axom::slam
{
/*!
 * \brief A static relation from \a FromSet to \a ToSet with per-element (variable) cardinality,
 *  reading offsets and indices from \c axom::Array buffers managed elsewhere.
 *
 * Matches \c make_variable_relation when the begins buffer uses \a FlatPosType.
 */
template <typename FromSet,
          typename ToSet,
          typename FlatPosType =
            detail::default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>>
using VariableRelation = StaticRelation<
  FlatPosType,
  typename ToSet::PositionType,
  policies::VariableCardinality<FlatPosType, policies::ArrayIndirection<FlatPosType, FlatPosType>>,
  policies::ArrayIndirection<FlatPosType, typename ToSet::PositionType>,
  FromSet,
  ToSet>;

/*!
 * \brief A static relation view from \a FromSet to \a ToSet with per-element
 *  (variable) cardinality, binding offsets and indices through \c axom::ArrayView.
 */
template <typename FromSet,
          typename ToSet,
          typename FlatPosType =
            detail::default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>>
using VariableRelationView = StaticRelation<
  FlatPosType,
  typename ToSet::PositionType,
  policies::VariableCardinality<FlatPosType, policies::ArrayViewIndirection<FlatPosType, FlatPosType>>,
  policies::ArrayViewIndirection<FlatPosType, typename ToSet::PositionType>,
  FromSet,
  ToSet>;

/*!
 * \brief A static relation from \a FromSet to \a ToSet with fixed cardinality \a N
 *  (each from-element maps to exactly N to-elements), reading indices from an \c axom::Array buffer managed elsewhere.
 */
template <typename FromSet,
          typename ToSet,
          int N,
          typename FlatPosType =
            detail::default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>>
using ConstantRelation =
  StaticRelation<FlatPosType,
                 typename ToSet::PositionType,
                 policies::ConstantCardinality<FlatPosType, policies::CompileTimeStride<FlatPosType, N>>,
                 policies::ArrayIndirection<FlatPosType, typename ToSet::PositionType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation view from \a FromSet to \a ToSet with fixed cardinality \a N,
 *  binding indices through \c axom::ArrayView.
 */
template <typename FromSet,
          typename ToSet,
          int N,
          typename FlatPosType =
            detail::default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>>
using ConstantRelationView =
  StaticRelation<FlatPosType,
                 typename ToSet::PositionType,
                 policies::ConstantCardinality<FlatPosType, policies::CompileTimeStride<FlatPosType, N>>,
                 policies::ArrayViewIndirection<FlatPosType, typename ToSet::PositionType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation from \a FromSet to \a ToSet with runtime constant cardinality,
 *  reading indices from an \c axom::Array buffer managed elsewhere.
 */
template <typename FromSet,
          typename ToSet,
          typename FlatPosType =
            detail::default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>>
using RuntimeConstantRelation =
  StaticRelation<FlatPosType,
                 typename ToSet::PositionType,
                 policies::ConstantCardinality<FlatPosType, policies::RuntimeStride<FlatPosType>>,
                 policies::ArrayIndirection<FlatPosType, typename ToSet::PositionType>,
                 FromSet,
                 ToSet>;

/*!
 * \brief A static relation view from \a FromSet to \a ToSet with runtime constant cardinality,
 *  binding indices through \c axom::ArrayView.
 */
template <typename FromSet,
          typename ToSet,
          typename FlatPosType =
            detail::default_flat_position_t<typename FromSet::PositionType, typename ToSet::PositionType>>
using RuntimeConstantRelationView =
  StaticRelation<FlatPosType,
                 typename ToSet::PositionType,
                 policies::ConstantCardinality<FlatPosType, policies::RuntimeStride<FlatPosType>>,
                 policies::ArrayViewIndirection<FlatPosType, typename ToSet::PositionType>,
                 FromSet,
                 ToSet>;

}  // namespace axom::slam
