// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef AXOM_MIR_EQUIZ_ALGORITHM_HPP_
#define AXOM_MIR_EQUIZ_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/mir/MIRAlgorithm.hpp"

#include <conduit/conduit.hpp>

namespace axom
{

namespace mir
{

class EquiZAlgorithm : public MIRAlgorithm
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

  EquiZAlgorithm() = default;
  virtual ~EquiZAlgorithm() = default;

  void setExecPolicy(RuntimePolicy policy) { m_execPolicy = policy; }

protected:
   /// Implement the EquiZ algorithm on a single domain.
   virtual void execute(const conduit::Node &topo,
                        const conduit::Node &coordset,
                        const conduit::Node &options,
                        conduit::Node &new_topo,
                        conduit::Node &new_coordset) override;

   
   /// Implement the EquiZ algorithm on a single domain for a given ExecSpace.
   template <typename ExecSpace>
   void executeImpl(const conduit::Node &topo,
                    const conduit::Node &coordset,
                    const conduit::Node &options,
                    conduit::Node &new_topo,
                    conduit::Node &new_coordset);

   RuntimePolicy m_execPolicy{RuntimePolicy::seq};
};

} // end namespace mir
} // end namespace axom

#endif
