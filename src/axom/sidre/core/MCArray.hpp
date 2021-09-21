// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SIDRE_MCARRAY_HPP_
#define SIDRE_MCARRAY_HPP_

#include "axom/core/Macros.hpp"  // for disable copy/assignment macro
#include "axom/core/utilities/Utilities.hpp"  // for memory allocation functions
#include "axom/core/MCArray.hpp"              // to inherit
#include "axom/core/Types.hpp"

#include "axom/slic/interface/slic.hpp"  // for slic logging macros

#include "View.hpp"    // for View definition
#include "Buffer.hpp"  // for Buffer definition
#include "Array.hpp"

// C/C++ includes
#include <cstring>  // for std::memcpy

namespace axom
{
namespace sidre
{
template <typename T>
using MCArray = Array<T, 2>;

} /* namespace sidre */
} /* namespace axom */

#endif /* SIDRE_MCARRAY_HPP_ */
