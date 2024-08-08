// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_CLIP_OPTIONS_HPP_
#define AXOM_MIR_CLIP_OPTIONS_HPP_

#include "axom/mir/Options.hpp"

namespace axom
{
namespace mir
{
namespace clipping
{

/**
 * \brief This class provides a kind of schema over the clipping options, as well
 *        as default values, and some utilities functions.
 */
template <typename ExecSpace>
class ClipOptions : public axom::mir::Options<ExecSpace>
{
public:
  /**
   * \brief Constructor
   *
   * \param nzones The total number of zones in the associated topology.
   * \param options The node that contains the clipping options.
   */
  ClipOptions(axom::IndexType nzones, const conduit::Node &options) : axom::mir::Options<ExecSpace>(nzones, options)
  { }

  /**
   * \brief Return the name of the field used for clipping.
   * \return The name of the field used for clipping.
   */
  std::string clipField() const
  {
    return options().fetch_existing("clipField").as_string();
  }

  /**
   * \brief Return the clip value.
   * \return The clip value.
   */
  float clipValue() const
  {
    return options().has_child("clipValue")
      ? options().fetch_existing("clipValue").to_float()
      : 0.f;
  }

  /**
   * \brief Return the name of the new color field to be created.
   * \return The name of the new color field to be created.
   */
  std::string colorField() const
  {
    std::string name("color");
    if(options().has_child("colorField"))
      name = options().fetch_existing("colorField").as_string();
    return name;
  }

  /**
   * \brief Whether the "inside" of the clipping field is selected.
   * \return 1 of the inside clipping is selected, false otherwise.
   */
  bool inside() const
  {
    return options().has_path("inside")
      ? (options().fetch_existing("inside").to_int() > 0)
      : true;
  }

  /**
   * \brief Whether the "outside" of the clipping field is selected.
   * \return 1 of the outside clipping is selected, false otherwise.
   */
  bool outside() const
  {
    return options().has_path("outside")
      ? (options().fetch_existing("outside").to_int() > 0)
      : false;
  }
private:
  /// Access the base class' options.
  const conduit::Node &options() const { return this->m_options; }
};

}  // end namespace clipping
}  // end namespace mir
}  // end namespace axom

#endif
