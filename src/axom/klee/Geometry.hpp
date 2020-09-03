// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_GEOMETRY_HPP
#define AXOM_KLEE_GEOMETRY_HPP

#include <string>

namespace axom { namespace klee {

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

private:
    std::string m_format;
    std::string m_path;
};

}}

#endif //AXOM_KLEE_GEOMETRY_HPP
