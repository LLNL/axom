// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Shape.hpp"

#include <algorithm>
#include <stdexcept>

namespace axom { namespace klee {

namespace {

// TODO: is there already a standard way to do this in axom? Is there a TPL
// already in use that does this? If not, is there a standard place to put
// utility functions like this?
/**
 * Check whether a given container contains a value
 *
 * \tparam Container the type of the container
 * \param container the container to check for the value
 * \param value the value to look for
 * \return whether the container contains the value
 */
template <typename Container>
bool contains(Container const &container,
        typename Container::value_type const &value) {
    using std::begin;
    using std::end;
    auto endIter = end(container);
    return std::find(begin(container), endIter, value) != endIter;
}
}

void Shape::setMaterialsReplaced(
        std::vector<std::string> const &materialsReplaced) {
    if (!m_materialsNotReplaced.empty()) {
        throw std::logic_error("Can't set list of materials to replace "
                               "when materials to not replace have already "
                               "been set");
    }
    m_materialsReplaced = materialsReplaced;
}

void Shape::setMaterialsNotReplaced(
        std::vector<std::string> const &materialsNotReplaced) {
    if (!m_materialsReplaced.empty()) {
        throw std::logic_error("Can't set list of materials to not replace "
                               "when materials to replace have already "
                               "been set");
    }
    m_materialsNotReplaced = materialsNotReplaced;
}

bool Shape::replaces(std::string const &material) const {
    if (!m_materialsReplaced.empty()) {
        return contains(m_materialsReplaced, material);
    } else if (!m_materialsNotReplaced.empty()) {
        return !contains(m_materialsNotReplaced, material);
    }
    return true;
}

void Shape::setName(std::string name) {
    m_name = std::move(name);
}

void Shape::setMaterial(std::string material) {
    m_material = std::move(material);
}

void Shape::setGeometry(Geometry geometry) {
    m_geometry = std::move(geometry);
}

}}
