// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_SHAPE_SET_HPP
#define AXOM_KLEE_SHAPE_SET_HPP

#include <string>
#include <vector>

#include "axom/klee/Dimensions.hpp"
#include "axom/klee/Shape.hpp"

namespace axom
{
namespace klee
{
/**
 * A ShapeSet represents a document in the common shape format.
 */
class ShapeSet
{
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
  std::vector<Shape> const &getShapes() const { return m_shapes; }

  /**
   * Set the file path from which this ShapeSet was created. This must be
   * set for resolvePath() to work.
   *
   * \param path the ShapeSet's path
   */
  void setPath(const std::string &path);

  /**
   * Get the path of the file from which this ShapeSet was created.
   *
   * \return the path of the file. Can be empty.
   */
  const std::string &getPath() const { return m_path; }

  /**
   * Resolves a path relative to the path of this ShapeSet.
   *
   * \param filePath the path to resolve
   * \return if the given path is absolute, then the given path. Otherwise,
   * the path is interpreted as being relative to the directory containing
   * this ShapeSet, and that is is returned.
   */
  std::string resolvePath(const std::string &filePath) const;

  /**
   * Sets the dimensions for all shapes in the ShapeSet.
   *
   * \param dimensions the dimension for all the shapes
   * \note This function must be called before calling \a getDimensions()
   */
  void setDimensions(Dimensions dimensions);

  /**
   * Returns the dimension of the ShapeSet.
   *
   * \pre Only valid after \a setDimensions() has been called on this instance
   * \sa setDimensions()
   */
  Dimensions getDimensions() const;

private:
  std::vector<Shape> m_shapes;
  std::string m_path;
  bool m_dimensionsHaveBeenSet {false};
  Dimensions m_dimensions;
};

}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEE_SHAPE_SET_HPP
