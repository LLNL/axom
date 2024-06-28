// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_UTILITIES_HPP_
#define AXOM_MIR_DISPATCH_UTILITIES_HPP_

namespace axom
{
namespace mir
{
namespace views
{

template <typename ... Dimensions>
constexpr int encode_dimensions(Dimensions... dims)
{
  return (... | dims);
}

template <typename ... Dimensions>
constexpr int select_dimensions(Dimensions... dims)
{
  return encode_dimensions((1 << dims)...);
}

constexpr bool dimension_selected(int encoded_dims, int dim)
{
  return encoded_dims & (1 << dim);
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
