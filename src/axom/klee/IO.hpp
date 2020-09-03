// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_IO_HPP
#define AXOM_KLEE_IO_HPP

#include <istream>

#include "axom/klee/ShapeSet.hpp"

namespace axom { namespace klee {

/**
 * Read a ShapeSet from an input stream.
 *
 * \param stream the stream from which to read the ShapeSet
 * \return the ShapeSet read from the stream
 * \throws runtime_error if the input is invalid
 */
ShapeSet readShapeSet(std::istream &stream);

}}

#endif //AXOM_KLEE_IO_HPP
