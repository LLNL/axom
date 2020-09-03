// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_SHAPE_SET_HPP
#define AXOM_KLEE_SHAPE_SET_HPP

#include <vector>

#include "axom/klee/Shape.hpp"

namespace axom { namespace klee {

/**
 * A ShapeSet represents a document in the common shape format.
 */
class ShapeSet {
public:
    /**
     * Set the shapes in this set.
     *
     * \param shapes all the shapes in this set
     */
    void setShapes(std::vector<Shape> shapes);

    /**
     * Get the shapes in this set.
     *
     * \return the shapes in this set
     */
    std::vector<Shape> const& getShapes() const {
        return m_shapes;
    }

private:
    std::vector<Shape> m_shapes;
};

}}


#endif //AXOM_KLEE_SHAPE_SET_HPP
