// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_OPTIONS_HPP_
#define AXOM_MIR_OPTIONS_HPP_

#include "axom/core.hpp"

#include <conduit/conduit.hpp>
#include <string>

namespace axom
{
namespace mir
{
// IDEA: maybe use inlet for this stuff.

/**
 * \brief This class provides a kind of schema over options, as well
 *        as default values, and some utilities functions.
 */
class Options
{
public:
  /**
   * \brief Constructor
   *
   * \param nzones The total number of zones in the associated topology.
   * \param options The node that contains the clipping options.
   */
  Options(const conduit::Node &options) : m_options(options) { }

  /**
   * \brief Return the name of the topology to make in the output.
   * \param default_value The name to use if the option is not defined.
   * \return The name of the topology to make in the output.
   */
  std::string topologyName(const std::string &default_value = std::string()) const
  {
    std::string name(default_value);
    if(m_options.has_child("topologyName"))
      name = m_options.fetch_existing("topologyName").as_string();
    return name;
  }

  /**
   * \brief Return the name of the coordset to make in the output.
   * \param default_value The name to use if the option is not defined.
   * \return The name of the coordset to make in the output.
   */
  std::string coordsetName(const std::string &default_value = std::string()) const
  {
    std::string name(default_value);
    if(m_options.has_child("coordsetName"))
      name = m_options.fetch_existing("coordsetName").as_string();
    return name;
  }

  /**
   * \brief Extract the names of the fields to process (and their output names) from the
   *        options or \a n_fields if the options do not contain fields.
   *
   * \param[out] f A map of the fields that will be processed, as well as their output name in the new fields.
   * \return True if the fields were present in the options. False otherwise.
   */
  bool fields(std::map<std::string, std::string> &f) const
  {
    bool retval = m_options.has_child("fields");
    f.clear();
    if(retval)
    {
      const conduit::Node &n_opt_fields = m_options.fetch_existing("fields");
      for(conduit::index_t i = 0; i < n_opt_fields.number_of_children(); i++)
      {
        if(n_opt_fields[i].dtype().is_string())
          f[n_opt_fields[i].name()] = n_opt_fields[i].as_string();
        else
          f[n_opt_fields[i].name()] = n_opt_fields[i].name();
      }
    }
    return retval;
  }

protected:
  const conduit::Node &m_options;  // A reference to the options node.
};
}  // end namespace mir
}  // end namespace axom

#endif
