// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/multimat/multimat.hpp"

#include <type_traits>

namespace axom
{
namespace multimat
{
namespace detail
{
// Select the parent's row view without stripping constness from the parent object.
template <typename Field>
using MMSubMapType =
  std::remove_const_t<std::conditional_t<std::is_const_v<Field>,
                                         typename Field::BiVarMapType::ConstSubMapType,
                                         typename Field::BiVarMapType::SubMapType>>;
}  // namespace detail

/**
 * Class for MultiMat 2D SubFields
 */

template <typename Field2DType>
class MMSubField2D : public detail::MMSubMapType<Field2DType>
{
public:
  using SubMapType = detail::MMSubMapType<Field2DType>;
  using SubSetType = typename SubMapType::IndexSetType;
  using SuperMapType = typename SubMapType::ParentMapType;
  using PositionType = typename SuperMapType::PositionType;
  using FirstPositionType = typename SuperMapType::FirstPositionType;
  using SecondPositionType = typename SuperMapType::SecondPositionType;
  using BiVarSetType = typename Field2DType::BiVarSetType;

  // Default Constructor
  MMSubField2D() = default;

  // Constructor
  AXOM_SUPPRESS_HD_WARN
  AXOM_HOST_DEVICE MMSubField2D(Field2DType* superfield,
                                FirstPositionType firstIndex,
                                bool indicesHaveIndirection = true)
    : SubMapType(static_cast<SuperMapType*>(superfield),
                 superfield->set()->elementRangeSet(firstIndex),
                 indicesHaveIndirection)
    , firstSetIndex(firstIndex)
  { }

  AXOM_HOST_DEVICE FirstPositionType getOuterIndex() const { return firstSetIndex; }

private:
  FirstPositionType firstSetIndex {-1};

};  //end class MMSubField2D

template <typename Field2DType, DataLayout DataLayoutT>
class MMSubField2DWrap : public MMSubField2D<Field2DType>
{ };

// specialization for Cell Dom
template <typename Field2DType>
class MMSubField2DWrap<Field2DType, DataLayout::CELL_DOM> : public MMSubField2D<Field2DType>
{
public:
  using SFB = MMSubField2D<Field2DType>;
  using typename SFB::FirstPositionType;
  using typename SFB::PositionType;
  using typename SFB::SecondPositionType;
  MMSubField2DWrap(Field2DType* superfield, FirstPositionType firstIndex, bool indirection = true)
    : SFB(superfield, firstIndex, indirection)
  { }

  DataLayout getDataLayout() { return DataLayout::CELL_DOM; }
  FirstPositionType cellId() const { return this->getOuterIndex(); }
  SecondPositionType matId(PositionType i) const { return this->index(i).second; }
};

// specialization for Mat Dom
template <typename Field2DType>
class MMSubField2DWrap<Field2DType, DataLayout::MAT_DOM> : public MMSubField2D<Field2DType>
{
public:
  using SFB = MMSubField2D<Field2DType>;
  using typename SFB::FirstPositionType;
  using typename SFB::PositionType;
  using typename SFB::SecondPositionType;
  MMSubField2DWrap(Field2DType* superfield, FirstPositionType firstIndex, bool indirection = true)
    : SFB(superfield, firstIndex, indirection)
  { }

  DataLayout getDataLayout() { return DataLayout::MAT_DOM; }
  FirstPositionType matId() const { return this->getOuterIndex(); }
  SecondPositionType cellId(PositionType i) const { return this->index(i).second; }
};

}  //end namespace multimat
}  //end namespace axom

///////////////////////////////////////////////////////////////////////////
