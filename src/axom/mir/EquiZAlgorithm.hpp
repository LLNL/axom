// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_EQUIZ_ALGORITHM_HPP_
#define AXOM_MIR_EQUIZ_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir/MIRAlgorithm.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace mir
{
/**
 * \accelerated
 * \brief Implements Meredith's Equi-Z algorithm on the GPU using Blueprint inputs/outputs.
 */
class EquiZAlgorithm : public MIRAlgorithm
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

  EquiZAlgorithm() = default;
  virtual ~EquiZAlgorithm() = default;

  /**
   * \brief Set the desired execution policy.
   *
   * \param policy The execution policy.
   *
   * \note The input Conduit nodes must have data located in a memory space
   *       that is compatible with the execution policy.
   */
  void setExecPolicy(RuntimePolicy policy) { m_execPolicy = policy; }

protected:
  /// Implement the EquiZ algorithm on a single domain.
  virtual void execute(const conduit::Node &topo,
                       const conduit::Node &coordset,
                       const conduit::Node &matset,
                       const conduit::Node &options,
                       conduit::Node &new_topo,
                       conduit::Node &new_coordset,
                       conduit::Node &new_matset) override;

  /// Implement the EquiZ algorithm on a single domain for a given ExecSpace.
  template <typename ExecSpace>
  void executeImpl(const conduit::Node &topo,
                   const conduit::Node &coordset,
                   const conduit::Node &matset,
                   const conduit::Node &options,
                   conduit::Node &new_topo,
                   conduit::Node &new_coordset,
                   conduit::Node &new_matset);

  RuntimePolicy m_execPolicy {RuntimePolicy::seq};
};

}  // end namespace mir
}  // end namespace axom

#endif
