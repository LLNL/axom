// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"

#include <utility>

namespace axom { namespace klee {

void Geometry::setFormat(std::string format) {
    m_format = std::move(format);
}

void Geometry::setPath(std::string path) {
    m_path = std::move(path);
}

}}
