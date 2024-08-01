// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_ALGORITHM_HPP_
#define AXOM_MIR_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/ArrayView.hpp"

#include <conduit/conduit.hpp>

#include <vector>
#include <string>

namespace axom
{
namespace mir
{
/**
 \brief Base class for Material Interface Reconstruction (MIR) algorithms.
 */
class MIRAlgorithm
{
public:
  MIRAlgorithm() = default;
  virtual ~MIRAlgorithm() = default;

  /**
    \brief Perform material interface reconstruction on the meshes supplied in the
           root node. Root can either be a mesh domain or a node that contains multiple
           domains.

    \param[in] root The root node that contains either a mesh or list of mesh
                    domains that contain a topology and matset to be used for MIR.
    \param[in] options A node that contains options that help govern MIR execution.

options:
  topology: main
  matset: matset
  new_topology: mirtopo
  new_coordset: mircoords
  new_matset: cleanmat
  fields:
    - temperature
    - pressure
  zones: [0,1,6,9]
  mapping: 0

    The "topology" option specifies which topology we'll reconstruct. It must have an associated matset.
    "new_topology" is the name of the topology that will be created in the output node.
    "new_coordset" is the name of the new coordset that will be created in the output node. If it is not provided then the name of the topology's coordset will be used.  
    "new_matset" is the name of the new matset that will be created in the output node. If it is not provided then the name of the topology's matset will be used.  
    "fields" is the name of the fields to map to the new topology. If fields is specified but empty, no fields will be mapped. If fields is not present then all fields will be mapped.
    "zones" is a list of zone indices from the topology that need to be reconstructed. If not present then all zones will be considered.
    "mapping" indicates whether we should include an original_element_numbers field on the new topology to indicate where each new zone came from in the original topology.

    \param[out] output A node that will contain the new entities.

    */
  virtual void execute(const conduit::Node &root,
                       const conduit::Node &options,
                       conduit::Node &output);

protected:
  /**
    \brief Perform material interface reconstruction on a single domain. Derived classes
           must implement this method and any device-specific coding gets handled under it.

    \param[in] topo The Conduit node containing the topology that will be used for MIR.
    \param[in] coordset The Conduit node containing the topology's coordset.
    \param[in] options The Conduit node containing the options that help govern MIR execution.

    \param[out] new_topo A Conduit node that will contain the new topology.
    \param[out] new_coordset A Conduit node that will contain the new coordset.
    
    */
  virtual void execute(const conduit::Node &topo,
                       const conduit::Node &coordset,
                       const conduit::Node &matset,
                       const conduit::Node &options,
                       conduit::Node &new_topo,
                       conduit::Node &new_coordset,
                       conduit::Node &new_matset) = 0;

  // Utility methods for derived types.
  void copyState(const conduit::Node &mesh, conduit::Node &destMesh) const;
  std::string topologyName(const conduit::Node &mesh,
                           const conduit::Node &options) const;
  std::string newTopologyName(const conduit::Node &mesh,
                              const conduit::Node &options) const;

  std::string newCoordsetName(const conduit::Node &mesh,
                              const conduit::Node &options) const;

  std::string matsetName(const conduit::Node &mesh,
                         const conduit::Node &options) const;
  std::string newMatsetName(const conduit::Node &mesh,
                            const conduit::Node &options) const;

  std::vector<std::string> fieldNames(const conduit::Node &mesh,
                                      const conduit::Node &options) const;

  const conduit::Node &topology(const conduit::Node &input,
                                const conduit::Node &options) const;
  const conduit::Node &matset(const conduit::Node &input,
                              const conduit::Node &options) const;

  // TODO: method for mapping element field to new topo
  // TODO: method for mapping vertex field to new topo
};

#if 0
class ElviraMIRAlgorithm : public MIRAlgorithm
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

  ElviraMIRAlgorithm() = default;
  virtual ~ElviraMIRAlgorithm() = default;

  void setExecPolicy(RuntimePolicy policy) { m_execPolicy = policy; }

protected:
   /// Implement the Elvira MIR algorithm on a single domain.
   virtual void execute(const conduit::Node &topo,
                        const conduit::Node &coordset,
                        const conduit::Node &options,
                        conduit::Node &new_topo,
                        conduit::Node &new_coordset) override;

   
   /// Implement the Elvira MIR algorithm on a single domain for a given ExecSpace.
   template <typename ExecSpace>
   void executeImpl(const conduit::Node &topo,
                    const conduit::Node &coordset,
                    const conduit::Node &options,
                    conduit::Node &new_topo,
                    conduit::Node &new_coordset);

   RuntimePolicy m_execPolicy{RuntimePolicy::seq};
};
#endif

}  // end namespace mir
}  // end namespace axom

#endif
