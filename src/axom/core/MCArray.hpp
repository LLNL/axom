// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MCARRAY_HPP_
#define AXOM_MCARRAY_HPP_

#include "axom/core/Array.hpp"  // for Array definition

namespace axom
{
template <typename T>
using MCArray = Array<T, 2>;

} /* namespace axom */

#endif /* AXOM_MCARRAY_HPP_ */
