// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_ELVIRA_ALGORITHM_HPP_
#define AXOM_MIR_ELVIRA_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

// Include these directly for now.
#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/utilities/ExtractZones.hpp"
#include "axom/mir/utilities/MakeZoneCenters.hpp"
#include "axom/mir/utilities/MergeMeshes.hpp"
#include "axom/mir/utilities/NodeToZoneRelationBuilder.hpp"
#include "axom/mir/utilities/RecenterField.hpp"
#include "axom/mir/utilities/SelectedZones.hpp"
#include "axom/mir/utilities/ZoneListBuilder.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/MaterialView.hpp"
#include "axom/mir/views/StructuredTopologyView.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <algorithm>
#include <string>

// Uncomment to save inputs and outputs.
#define AXOM_ELVIRA_DEBUG

#if defined(AXOM_ELVIRA_DEBUG)
  #include <conduit/conduit_relay_io.hpp>
  #include <conduit/conduit_relay_io_blueprint.hpp>
#endif

namespace axom
{
namespace mir
{
namespace elvira
{

/*!
 * \brief Compute the jacobian from an input stencil of x,y,z coordinates derived
 *        from zone centroids.
 *
 * \param xcst A stencil of X coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param ycst A stencil of Y coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param zcst A stencil of Z coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param ndims The number of dimensions (2 or 3).
 * \param[out] jac A 3x3 Jacobian matrix.
 *
 * \note  Adapted from J. Grandy's BasicStencil.cc:52 
 */
AXOM_HOST_DEVICE
void computeJacobian(const double *xcst,
                     const double *ycst,
                     const double *zcst,
                     int ndims,
                     double jac[3][3]);

/*!
 * \brief Multiply the normal times the Jacobian and normalize.
 *
 * \param normal The normal to transform, nx,ny,nz
 * \param jac    The jacobian to multiply by
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
void transform(double *normal, const double jac[3][3]);

/*!
 * \brief Get the size of the stencil based on dimension.
 *
 * \param ndims The stencil size.
 *
 * \return The stencil size for dimension \a ndims.
 */
AXOM_HOST_DEVICE
constexpr int getStencilSize(int ndims)
{
  return (ndims == 3) ? 27 : ((ndims == 2) ? 9 : 3);
}

/*! 
 * \brief Base template for a class that invokes ELVIRA on various dimension data.
 */
template <int NDIMS>
struct elvira {};

/*!
 * \brief 2D specialization that calls elvira2xy to make normals.
 */
template <>
struct elvira<2>
{
  static constexpr int NDIMS = 2;

  /*!
   * \brief Create normals for the material interface fragments in the selected zone.
   *
   * \param matCount The number of materials in the current zone.
   * \param fragmentVFStencilStart The start of the material volume fraction stencil data for all fragments in a zone.
   * \param fragmentVectorsStart The start of the normals for all fragments in a zone.
   * \param iskip The material index to skip.
   *
   * \note Calling this function will update some vectors in the \a fragmentVectorsStart.
   */
  AXOM_HOST_DEVICE
  static void execute(int matCount,
                      const double *fragmentVFStencilStart,
                      double *fragmentVectorsStart,
                      int iskip);
};

/*!
 * \brief 3D specialization that calls elvira3d to make normals.
 */
template <>
struct elvira<3>
{
  static constexpr int NDIMS = 3;

  /*!
   * \brief Create normals for the material interface fragments in the selected zone.
   *
   * \param matCount The number of materials in the current zone.
   * \param fragmentVFStencilStart The start of the material volume fraction stencil data for all fragments in a zone.
   * \param fragmentVectorsStart The start of the normals for all fragments in a zone.
   * \param iskip The material index to skip.
   *
   * \note Calling this function will update some vectors in the \a fragmentVectorsStart.
   */
  AXOM_HOST_DEVICE
  static void execute(int matCount,
                      const double *fragmentVFStencilStart,
                      double *fragmentVectorsStart,
                      int iskip);
};

}  // namespace elvira
//------------------------------------------------------------------------------

/*!
 * \brief Implements Elvira algorithm for structured meshes.
 */
template <typename ExecSpace, typename IndexPolicy, typename CoordsetView, typename MatsetView>
class ElviraAlgorithm : public axom::mir::MIRAlgorithm
{
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

public:
  using TopologyView = axom::mir::views::StructuredTopologyView<IndexPolicy>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  ElviraAlgorithm(const TopologyView &topoView,
                  const CoordsetView &coordsetView,
                  const MatsetView &matsetView)
    : axom::mir::MIRAlgorithm()
    , m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~ElviraAlgorithm() = default;

  /*!
   * \brief Perform material interface reconstruction on a single domain.
   *
   * \param[in] n_topo The Conduit node containing the topology that will be used for MIR.
   * \param[in] n_coordset The Conduit node containing the coordset.
   * \param[in] n_fields The Conduit node containing the fields.
   * \param[in] n_matset The Conduit node containing the matset.
   * \param[in] n_options The Conduit node containing the options that help govern MIR execution.
   *
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   * \param[out] n_newMatset A Conduit node that will contain the new matset.
   * 
   */
  virtual void executeDomain(const conduit::Node &n_topo,
                             const conduit::Node &n_coordset,
                             const conduit::Node &n_fields,
                             const conduit::Node &n_matset,
                             const conduit::Node &n_options,
                             conduit::Node &n_newTopo,
                             conduit::Node &n_newCoordset,
                             conduit::Node &n_newFields,
                             conduit::Node &n_newMatset) override
  {
    namespace bputils = axom::mir::utilities::blueprint;

    AXOM_ANNOTATE_SCOPE("ElviraAlgorithm");

    // Copy the options to make sure they are in the right memory space.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

    // Get selected zones from the options.
    bputils::SelectedZones<ExecSpace> selectedZones(
      m_topologyView.numberOfZones(),
      n_options_copy);
    const auto selectedZonesView = selectedZones.view();

    // Partition the selected zones into clean, mixed lists.
    axom::Array<axom::IndexType> cleanZones, mixedZones;
    bputils::ZoneListBuilder<ExecSpace, TopologyView, MatsetView> zlb(
      m_topologyView,
      m_matsetView);
    zlb.execute(selectedZonesView, cleanZones, mixedZones);

#if 1
    conduit::Node n_result;
    n_result["selectedZonesView"].set(selectedZonesView.data(), selectedZonesView.size());
    n_result["cleanZones"].set(cleanZones.data(), cleanZones.size());
    n_result["mixedZones"].set(mixedZones.data(), mixedZones.size());
    conduit::relay::io::save(n_result, "zonelist.yaml", "yaml");
#endif

    if(cleanZones.size() > 0 && mixedZones.size() > 0)
    {
      // Some clean, some mixed.

      // TODO: extract clean zones.

      processMixedZones(mixedZones.view(),
                        n_topo,
                        n_coordset,
                        n_fields,
                        n_matset,
                        n_options,
                        n_newTopo,
                        n_newCoordset,
                        n_newFields,
                        n_newMatset);

      // Merge meshes
    }
    else if(mixedZones.size() > 0)
    {
      // All mixed.
      processMixedZones(mixedZones.view(),
                        n_topo,
                        n_coordset,
                        n_fields,
                        n_matset,
                        n_options,
                        n_newTopo,
                        n_newCoordset,
                        n_newFields,
                        n_newMatset);
    }
    else if(cleanZones.size() > 0)
    {
      // All clean. We just pass the inputs through.
      n_newTopo.set_external(n_topo);
      n_newCoordset.set_external(n_coordset);
      n_newFields.set_external(n_fields);
      n_newMatset.set_external(n_matset);
    }
  }

  void processMixedZones(const axom::ArrayView<axom::IndexType> mixedZonesView,
                         const conduit::Node &n_topo,
                         const conduit::Node &n_coordset,
                         const conduit::Node &n_fields,
                         const conduit::Node &n_matset,
                         const conduit::Node &n_options,
                         conduit::Node &n_newTopo,
                         conduit::Node &n_newCoordset,
                         conduit::Node &n_newFields,
                         conduit::Node &n_newMatset)
  {
    AXOM_ANNOTATE_SCOPE("processMixedZones");
    namespace bputils = axom::mir::utilities::blueprint;
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    constexpr int ndims = TopologyView::dimension();

#ifdef AXOM_ELVIRA_DEBUG
    // Let's output the normals
    conduit::Node *n_result = new conduit::Node;
#endif

    // Count the number of fragments we'll make for the mixed zones.
    // We also
    AXOM_ANNOTATE_BEGIN("counting");
    RAJA::ReduceSum<reduce_policy, axom::IndexType> num_reduce(0);

    const auto nzones = mixedZonesView.size();
    axom::Array<axom::IndexType> matCount(nzones, nzones, allocatorID);
    axom::Array<axom::IndexType> matZone(nzones, nzones, allocatorID);
    auto matCountView = matCount.view();
    auto matZoneView = matZone.view();

    const auto matsetView = m_matsetView;
    axom::for_all<ExecSpace>(
      mixedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = mixedZonesView[szIndex];
        const auto nmats = matsetView.numberOfMaterials(zoneIndex);
        matCountView[szIndex] = nmats;
        matZoneView[szIndex] = zoneIndex;
        num_reduce += nmats;
      });
    const auto numFragments = num_reduce.get();
    AXOM_ANNOTATE_END("counting");

#ifdef AXOM_ELVIRA_DEBUG
    conduit::Node &n_group1 = n_result->operator[]("group1");
    n_group1["mixedZones"].set(mixedZonesView.data(), mixedZonesView.size());
    n_group1["matZone"].set(matZoneView.data(), matZoneView.size());
    n_group1["matCount"].set(matCountView.data(), matCountView.size());
    n_group1["numFragments"].set(numFragments);
#endif
    //--------------------------------------------------------------------------
    // Sort the zones by the mat count. This should make adjacent zones in the
    // list more likely to have the same number of materials.
    //
    // NOTE: When I look at zone numbers in ranges where they have been sorted
    //       by matCount, the zone numbers are wildly out of order. Do we need
    //       a different sort? Or, to sort the subranges of zones?
    AXOM_ANNOTATE_BEGIN("sorting");
    RAJA::sort_pairs<loop_policy>(RAJA::make_span(matCountView.data(), nzones),
                                  RAJA::make_span(matZoneView.data(), nzones));
    AXOM_ANNOTATE_END("sorting");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("offsets");
    axom::Array<axom::IndexType> matOffset(nzones, nzones, allocatorID);
    auto matOffsetView = matOffset.view();
    axom::exclusive_scan<ExecSpace>(matCountView, matOffsetView);
    AXOM_ANNOTATE_END("offsets");

#ifdef AXOM_ELVIRA_DEBUG
    conduit::Node &n_group2 = n_result->operator[]("group2");
    n_group2["matZone"].set(matZoneView.data(), matZoneView.size());
    n_group2["matCount"].set(matCountView.data(), matCountView.size());
    n_group2["offsets"].set(matOffsetView.data(), matOffsetView.size());
#endif

    //--------------------------------------------------------------------------
    // NOTE: this makes zone centers for all zones.
    AXOM_ANNOTATE_BEGIN("centroids");
    bputils::MakeZoneCenters<ExecSpace, TopologyView, CoordsetView> zc(
      m_topologyView,
      m_coordsetView);
    conduit::Node n_zcfield;
    zc.execute(n_topo, n_coordset, n_zcfield);
    using zc_value = typename CoordsetView::value_type;
    axom::ArrayView<zc_value> xview, yview, zview;
    xview = bputils::make_array_view<zc_value>(n_zcfield["values/x"]);
    yview = bputils::make_array_view<zc_value>(n_zcfield["values/y"]);
    if(n_zcfield.has_path("values/z"))
    {
      zview = bputils::make_array_view<zc_value>(n_zcfield["values/z"]);
    }
    AXOM_ANNOTATE_END("centroids");
#ifdef AXOM_ELVIRA_DEBUG
    conduit::Node &n_group3 = n_result->operator[]("group3");
    n_group3["xview"].set(xview.data(), xview.size());
    n_group3["yview"].set(yview.data(), yview.size());
    n_group3["zview"].set(zview.data(), zview.size());
#endif

    //--------------------------------------------------------------------------
    // Retrieve stencil data for each zone material.
    AXOM_ANNOTATE_BEGIN("stencil");
    constexpr int StencilSize = elvira::getStencilSize(ndims);
    const auto numFragmentsStencil = numFragments * StencilSize;

    axom::Array<double> fragmentVFStencil(numFragmentsStencil,
                                          numFragmentsStencil,
                                          allocatorID);
    auto fragmentVFStencilView = fragmentVFStencil.view();

    // Sorted material ids for each zone.
    axom::Array<typename MatsetView::IndexType> sortedMaterialIds(numFragments,
                                                                  numFragments,
                                                                  allocatorID);
    auto sortedMaterialIdsView = sortedMaterialIds.view();

    // Coordinate stencil data for each zone.
    const auto nzonesStencil = nzones * StencilSize;
    axom::Array<double> xcStencil(nzonesStencil, nzonesStencil, allocatorID);
    axom::Array<double> ycStencil(nzonesStencil, nzonesStencil, allocatorID);
    axom::Array<double> zcStencil(nzonesStencil, nzonesStencil, allocatorID);
    auto xcStencilView = xcStencil.view();
    auto ycStencilView = ycStencil.view();
    auto zcStencilView = zcStencil.view();

    // Traverse the selected zones based on how many materials there are in a zone.
    TopologyView deviceTopologyView(m_topologyView);
    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        const auto matCount = matCountView[szIndex];
        // Where to begin writing this zone's fragment data.
        const auto offset = matOffsetView[szIndex];

        // Get materials for this zone from the matset.
        typename MatsetView::IDList ids;
        typename MatsetView::VFList vfs;
        matsetView.zoneMaterials(zoneIndex, ids, vfs);

        // Sort the materials by the volume fraction. Save sorted ids in sortedMaterialIdsView.
        axom::utilities::sort_multiple(vfs.data(), ids.data());
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          const auto fragmentIndex = offset + m;
          sortedMaterialIdsView[fragmentIndex] = ids[m];
        }

        // Retrieve the stencil data from neighbor zones.
        auto logical = deviceTopologyView.indexing().IndexToLogicalIndex(zoneIndex);
        for(int si = 0; si < StencilSize; si++)
        {
          // Stencil neighbor logical index.
          typename TopologyView::LogicalIndex neighbor(logical);
          // Neighbor offsets are in (-1, 0, 1) that get added to current zone's logical coordinate.
          const int neighborOffset[3] = {(si % 3) - 1, ((si % 9) / 3) - 1, (si / 9) - 1};
          for(int d = 0; d < ndims; d++)
          {
            neighbor[d] += neighborOffset[d];
          }

          const auto coordIndex = szIndex * StencilSize + si;

          if(deviceTopologyView.indexing().contains(neighbor))
          {
            const auto neighborIndex = static_cast<typename MatsetView::ZoneIndex>(
              deviceTopologyView.indexing().LogicalIndexToIndex(neighbor));

            for(axom::IndexType m = 0; m < matCount; m++)
            {
              // Ask the neighbor zone for how much of the current material it contains.
              const auto fragmentIndex = offset + m;
              typename MatsetView::FloatType vf = 0;
              matsetView.zoneContainsMaterial(
                neighborIndex,
                sortedMaterialIdsView[fragmentIndex],
                vf);

              // Store the vf into the stencil for the current material.
              const auto destIndex = fragmentIndex * StencilSize + si;
              fragmentVFStencilView[destIndex] = static_cast<double>(vf);
            }

            // coord stencil
            xcStencilView[coordIndex] = xview[neighborIndex];
            ycStencilView[coordIndex] = yview[neighborIndex];
            zcStencilView[coordIndex] = zview.empty() ? 0. : zview[neighborIndex];
          }
          else
          {
            // All of the material contributions for the neighbor are 0.
            for(axom::IndexType m = 0; m < matCount; m++)
            {
              // Store the vf into the stencil for the current material.
              const auto fragmentIndex = offset + m;
              const auto destIndex = fragmentIndex * StencilSize + si;
              fragmentVFStencilView[destIndex] = 0.;
            }

            // coord stencil
            xcStencilView[coordIndex] = 0.;
            ycStencilView[coordIndex] = 0.;
            zcStencilView[coordIndex] = 0.;
          }
        }
      });
    AXOM_ANNOTATE_END("stencil");

    AXOM_ANNOTATE_BEGIN("vectors");
    constexpr int numVectorComponents = 3;
    const auto vecSize = numFragments * numVectorComponents;
    axom::Array<double> fragmentVectors(vecSize, vecSize, allocatorID);
    auto fragmentVectorsView = fragmentVectors.view();

#if defined(AXOM_ELVIRA_DEBUG) && !defined(AXOM_DEVICE_CODE)
    conduit::Node *n_group4 = &(n_result->operator[]("group4"));
#endif

    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        const auto matCount = matCountView[szIndex];
        // Where to begin writing this zone's fragment data.
        const auto offset = matOffsetView[szIndex];

        // The following lines are adapted from mira1.c:262

        // Compute Jacobian here since we have coordinate stencil data.
        double jac[3][3];
        const auto coordIndex = szIndex * StencilSize;
        const double *xcStencil = xcStencilView.data() + coordIndex;
        const double *ycStencil = ycStencilView.data() + coordIndex;
        const double *zcStencil = zcStencilView.data() + coordIndex;
        elvira::computeJacobian(xcStencil, ycStencil, zcStencil, ndims, jac);

        // The starting addresses for fragments in the current zone.
        const double *fragmentVFStencilStart = fragmentVFStencilView.data() + offset * StencilSize;
        double *fragmentVectorsStart = fragmentVectorsView.data() + offset * numVectorComponents;

        // Produce normal for each material in this zone.
        int iskip = 0;
        elvira::elvira<ndims>::execute(
               matCount,
               fragmentVFStencilStart,
               fragmentVectorsStart,
               iskip);

#if defined(AXOM_ELVIRA_DEBUG) && !defined(AXOM_DEVICE_CODE)
        {
        conduit::Node &n_thisZone = n_group4->append();
        n_thisZone["szIndex"] = szIndex;
        n_thisZone["zone"] = zoneIndex;
        n_thisZone["matCount"] = matCount;
        n_thisZone["offset"] = offset;
        n_thisZone["xcStencil"].set(xcStencil, StencilSize);
        n_thisZone["ycStencil"].set(ycStencil, StencilSize);
        n_thisZone["zcStencil"].set(zcStencil, StencilSize);
        n_thisZone["jacobian"].set(&jac[0][0], 9);
        conduit::Node &n_mats = n_thisZone["mats"];
        const double *vf = fragmentVFStencilStart;
        double *n = fragmentVectorsStart;
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          conduit::Node &n_thismat = n_mats.append();
          n_thismat["mat"] = sortedMaterialIdsView[offset + m];
          n_thismat["stencil"].set(vf, StencilSize);
          n_thismat["normal"].set(n, 3);
          vf += StencilSize;
          n += numVectorComponents;
        }
        }
#endif

        // Transform the normals.
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          double *normal =
            fragmentVectorsView.data() + ((offset + m) * numVectorComponents);
          elvira::transform(normal, jac);
        }
      });
    AXOM_ANNOTATE_END("vectors");

#ifdef AXOM_ELVIRA_DEBUG
    conduit::relay::io::save(*n_result, "elvira.yaml", "yaml");
    delete n_result;
#endif

    // intermediate...

    // Now that we have all of the normals for each zone, clip the zones using those normals, making a fragment
    // shape that has the right VF. then save the shape, take the remaining piece and apply the next normal to it.

    // Add in clean zones - do we end up with a vector of polyhedral zones here that might be easier to pass into
    // the mapping stage? Or, do we make a BG topology for the polyhedra and make sure the mapper can eat that?

    // Make matset for combined mesh
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif
