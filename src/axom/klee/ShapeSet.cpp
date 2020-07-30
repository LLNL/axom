// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/ShapeSet.hpp"

#include <utility>

#include "conduit.hpp"

namespace axom { namespace klee {

void ShapeSet::setShapes(std::vector<Shape> shapes) {
    m_shapes = std::move(shapes);
}

}}
