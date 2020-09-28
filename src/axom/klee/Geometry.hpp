// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_GEOMETRY_HPP
#define AXOM_KLEE_GEOMETRY_HPP

#include <memory>
#include <string>

namespace axom { namespace klee {

class GeometryOperator;

/**
 * Represents the geometry specified in a Shape.
 */
class Geometry {
public:
    /**
     * Get the format in which the geometry is specified.
     *
     * \return the format of the shape
     */
    const std::string &getFormat() const {
        return m_format;
    }

    /**
     * Set the format.
     *
     * \param format the geometry's format
     */
    void setFormat(std::string format);

    /**
     * Get the path at which to find the specification of the geometry
     *
     * \return the path to the geometry file
     */
    const std::string &getPath() const {
        return m_path;
    }

    /**
     * Set the path to the geometry file.
     *
     *  \param path the file's path
     */
    void setPath(std::string path);

    /**
     * Get a GeometryOperator to apply to this geometry. Can be null.
     *
     * \return a potentially null operator to apply to the geometry
     */
    std::shared_ptr<GeometryOperator const> const &getGeometryOperator() const {
        return m_operator;
    }

    /**
     * Set a GeometryOperator to apply to this geometry.
     *
     * \param op the operator to apply
     */
    void setGeometryOperator(std::shared_ptr<GeometryOperator const> const &op) {
        m_operator = op;
    }

    /**
     * Set the initial dimensions of this geometry, before any operators
     * have taken effect.
     *
     * \param initialDimensions the initial dimensions. Must be 2 or 3.
     */
    void setInitialDimensions(int initialDimensions) {
        m_initialDimensions = initialDimensions;
    }

    /**
     * Set the initial dimensions of this geometry
     *
     * \return the initial dimensions of this geometry
     */
    int getInitialDimensions() const {
        return m_initialDimensions;
    }

    /**
     * The dimensions of this geometry, after applying any operators.
     *
     * \return the geometry's dimension
     */
    int getDimensions() const;

private:
    int m_initialDimensions;
    std::string m_format;
    std::string m_path;
    std::shared_ptr<const GeometryOperator> m_operator;
};

}}

#endif //AXOM_KLEE_GEOMETRY_HPP
