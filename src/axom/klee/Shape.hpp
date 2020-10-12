// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_SHAPE_HPP
#define AXOM_KLEE_SHAPE_HPP

#include <string>
#include <vector>

#include "axom/klee/Geometry.hpp"

namespace axom
{
namespace klee
{
/**
 * A Shape describes a region of space, a material it is made of, and a list
 * of materials which it does/doesn't replace.
 */
class Shape
{
public:
  /**
   * Get the name of this shape.
   * \return the shape's name
   */
  const std::string &getName() const { return m_name; }

  /**
   * Set the name of the shape.
   * \param name the shape's name
   */
  void setName(std::string name);

  /**
   * Get the material this shape is made of.
   * \return the shape's material.
   */
  const std::string &getMaterial() const { return m_material; }

  /**
   * Set the material this shape is made of.
   * \return the shape's material.
   */
  void setMaterial(std::string material);

  /**
   * Set the list of materials which can be replaced. By default, all
   * materials are replaced.
   *
   * \param materialsReplaced the list of materials which this shape can
   * replace.
   * \throws std::logic_error if this shape has a list of materials that
   * cannot be replaced.
   */
  void setMaterialsReplaced(const std::vector<std::string> &materialsReplaced);

  /**
   * Set the list of materials which cannot be replaced. By default, all
   * materials are replaced.
   *
   * \param materialsNotReplaced the list of materials which this shape
   * cannot replace.
   * \throws std::logic_error if this shape has a list of materials that
   * can be replaced.
   */
  void setMaterialsNotReplaced(const std::vector<std::string> &materialsNotReplaced);

  /**
   * Check whether this shape can replace the given material (within the
   * space defined by this shape).
   *
   * \param material the material to check
   * \return whether this shape replaces the given material
   */
  bool replaces(const std::string &material) const;

  /**
   * Get the description fo the geometry for this shape.
   *
   * \return the shape's geometry
   */
  const Geometry &getGeometry() const { return m_geometry; }

  /**
   * Set the shape's geometry
   * \param geometry a description of the shape's geometry
   */
  void setGeometry(Geometry geometry);

private:
  std::string m_name;
  std::string m_material;
  std::vector<std::string> m_materialsReplaced;
  std::vector<std::string> m_materialsNotReplaced;
  Geometry m_geometry;
};

}  // namespace klee
}  // namespace axom

#endif
