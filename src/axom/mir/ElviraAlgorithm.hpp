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
#include "axom/mir/utilities/MakeZoneVolumes.hpp"
#include "axom/mir/utilities/MergeMeshes.hpp"
#include "axom/mir/utilities/PrimalAdaptor.hpp"
#include "axom/mir/utilities/NodeToZoneRelationBuilder.hpp"
#include "axom/mir/utilities/RecenterField.hpp"
#include "axom/mir/utilities/SelectedZones.hpp"
#include "axom/mir/utilities/ZoneListBuilder.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"
#include "axom/mir/utilities/utilities.hpp"
//#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/MaterialView.hpp"
#include "axom/mir/views/StructuredTopologyView.hpp"

#include "axom/primal/operators/compute_bounding_box.hpp"
#include "axom/primal/operators/clip.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Plane.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_relay_io.hpp>
#include <conduit/conduit_relay_io_blueprint.hpp>

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
struct elvira
{ };

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
 * \brief Base template for building a new Conduit mesh.
 */
template <typename ExecSpace,
          typename CoordsetView,
          typename TopologyView,
          typename MatsetView,
          typename PolygonShape,
          int NDIMS>
struct TopologyBuilder
{ };

/*!
 * \brief Template specialization for making a new 2D Conduit polygon mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam CoordsetView The view that wraps the coordset data.
 * \tparam TopologyView The view that wraps the topology data.
 * \tparam MatsetView The view that wraps the matset data.
 * \tparam PolygonShape An axom::primal::Polygon (with various template args filled in).
 */
template <typename ExecSpace,
          typename CoordsetView,
          typename TopologyView,
          typename MatsetView,
          typename PolygonShape>
class TopologyBuilder<ExecSpace, CoordsetView, TopologyView, MatsetView, PolygonShape, 2>
{
  // The way we build fragments from quads, we should not get any with more than 5 sides.
  static constexpr int MAX_POINTS_PER_FRAGMENT = 5;

  using CoordType = typename CoordsetView::value_type;
  using ConnectivityType = typename TopologyView::ConnectivityType;
  using MaterialID = typename MatsetView::MaterialID;
  using MaterialVF = typename MatsetView::FloatType;

public:
  /*!
   * \brief Create the Blueprint nodes for the output coordset, topology, fields, and matset.
   *
   * \note This is a host-only function. We allocate fixed size arrays for the coordset and
   *       connectivity that will have gaps that consumers will need to skip over using
   *       the provided sizes/offsets.
   */
  void allocate(axom::IndexType numFragments,
                conduit::Node &n_coordset,
                conduit::Node &n_topology,
                conduit::Node &n_fields,
                conduit::Node &n_matset)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    const auto numCoordValues = numFragments * MAX_POINTS_PER_FRAGMENT;

    // Set up coordset and allocate data arrays.
    // Note that we overallocate the number of nodes to numCoordValues.
    n_coordset["type"] = "explicit";
    n_coordset["values/x"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/x"].set(
      conduit::DataType(bputils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_x = bputils::make_array_view<CoordType>(n_coordset["values/x"]);
    n_coordset["values/y"].set_allocator(c2a.getConduitAllocatorID());
    n_coordset["values/y"].set(
      conduit::DataType(bputils::cpp2conduit<CoordType>::id, numCoordValues));
    m_view.m_y = bputils::make_array_view<CoordType>(n_coordset["values/y"]);

    axom::mir::utilities::fill<ExecSpace>(m_view.m_x, CoordType(0));
    axom::mir::utilities::fill<ExecSpace>(m_view.m_y, CoordType(0));

    // Set up connectivity and allocate data arrays.
    n_topology["type"] = "unstructured";
    n_topology["elements/shape"] = "polygonal";
    conduit::Node &n_conn = n_topology["elements/connectivity"];
    n_conn.set_allocator(c2a.getConduitAllocatorID());
    n_conn.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id,
                                 numFragments * MAX_POINTS_PER_FRAGMENT));
    m_view.m_connectivity = bputils::make_array_view<ConnectivityType>(n_conn);
    axom::mir::utilities::fill<ExecSpace>(m_view.m_connectivity,
                                          ConnectivityType(0));

    conduit::Node &n_sizes = n_topology["elements/sizes"];
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id,
                                  numFragments));
    m_view.m_sizes = bputils::make_array_view<ConnectivityType>(n_sizes);

    conduit::Node &n_offsets = n_topology["elements/offsets"];
    n_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id,
                                    numFragments));
    m_view.m_offsets = bputils::make_array_view<ConnectivityType>(n_offsets);

    // Make new fields.
    n_fields["originalZones/topology"] = n_topology.name();
    n_fields["originalZones/association"] = "element";
    conduit::Node &n_orig_zones = n_fields["originalZones/values"];
    n_orig_zones.set_allocator(c2a.getConduitAllocatorID());
    n_orig_zones.set(conduit::DataType(bputils::cpp2conduit<ConnectivityType>::id,
                                       numFragments));
    m_view.m_original_zones =
      bputils::make_array_view<ConnectivityType>(n_orig_zones);

    conduit::Node &n_normal = n_fields["normal"];
    n_normal["topology"] = n_topology.name();
    n_normal["association"] = "element";
    conduit::Node &n_x = n_normal["values/x"];
    conduit::Node &n_y = n_normal["values/y"];
    n_x.set_allocator(c2a.getConduitAllocatorID());
    n_x.set(conduit::DataType(bputils::cpp2conduit<double>::id, numFragments));
    m_view.m_norm_x = bputils::make_array_view<double>(n_x);
    n_y.set_allocator(c2a.getConduitAllocatorID());
    n_y.set(conduit::DataType(bputils::cpp2conduit<double>::id, numFragments));
    m_view.m_norm_y = bputils::make_array_view<double>(n_y);

    // Set up new matset. All of the sizes are numFragments because we're making clean zones.
    n_matset["topology"] = n_topology.name();
    conduit::Node &n_volume_fractions = n_matset["volume_fractions"];
    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_volume_fractions.set(
      conduit::DataType(bputils::cpp2conduit<MaterialVF>::id, numFragments));
    m_view.m_volume_fractions =
      bputils::make_array_view<MaterialVF>(n_volume_fractions);

    conduit::Node &n_material_ids = n_matset["material_ids"];
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set(
      conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_material_ids = bputils::make_array_view<MaterialID>(n_material_ids);

    conduit::Node &n_indices = n_matset["indices"];
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set(
      conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_indices = bputils::make_array_view<MaterialID>(n_indices);

    conduit::Node &n_mat_sizes = n_matset["sizes"];
    n_mat_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_mat_sizes.set(
      conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_sizes = bputils::make_array_view<MaterialID>(n_mat_sizes);
    axom::mir::utilities::fill<ExecSpace>(m_view.m_mat_sizes, MaterialID(0));

    conduit::Node &n_mat_offsets = n_matset["offsets"];
    n_mat_offsets.set_allocator(c2a.getConduitAllocatorID());
    n_mat_offsets.set(
      conduit::DataType(bputils::cpp2conduit<MaterialID>::id, numFragments));
    m_view.m_mat_offsets = bputils::make_array_view<MaterialID>(n_mat_offsets);
  }

  /*!
   * \brief A view that can be used on device to encode a shape into Blueprint arrays.
   */
  struct View
  {
    /*!
     * \brief Encode a polygon from MIR into the views that represent the Conduit data.
     *
     * \param zoneIndex The original zone that yielded the fragment.
     * \param fragmentOffset The offset into the per-fragment data.
     * \param shape The 2D primal polygon we're encoding.
     * \param matId The material id for the fragment.
     * \param normal The normal for the fragment.
     */
    AXOM_HOST_DEVICE
    void addShape(axom::IndexType zoneIndex,
                  axom::IndexType fragmentOffset,
                  const PolygonShape &shape,
                  int matId,
                  const double *normal) const
    {
      const int nverts = shape.numVertices();
      SLIC_ASSERT(nverts <= MAX_POINTS_PER_FRAGMENT);

      // Copy coordinates into coordinate arrays. We might end up with unreferenced coordinates for now.
      auto coordOffset = fragmentOffset * MAX_POINTS_PER_FRAGMENT;
      for(int i = 0; i < nverts; i++)
      {
        m_x[coordOffset + i] = shape[i][0];
        m_y[coordOffset + i] = shape[i][1];
      }

      // Make connectivity.
      auto connOffset = fragmentOffset * MAX_POINTS_PER_FRAGMENT;
      for(int i = 0; i < nverts; i++)
      {
        auto idx = static_cast<ConnectivityType>(connOffset + i);
        m_connectivity[idx] = idx;
      }
      m_sizes[fragmentOffset] = static_cast<ConnectivityType>(nverts);
      m_offsets[fragmentOffset] = connOffset;

      // Save fields.
      m_original_zones[fragmentOffset] = zoneIndex;
      m_norm_x[fragmentOffset] = normal[0];
      m_norm_y[fragmentOffset] = normal[1];

      // Save material data.
      m_volume_fractions[fragmentOffset] = static_cast<MaterialVF>(1.);
      m_material_ids[fragmentOffset] = static_cast<MaterialID>(matId);
      m_mat_sizes[fragmentOffset] = static_cast<MaterialID>(1);
      m_mat_offsets[fragmentOffset] = fragmentOffset;
      m_mat_indices[fragmentOffset] = fragmentOffset;
    }

    // These ArrayView objects expose data that we've allocated in Conduit nodes to represent the new mesh.

    // New coordset data.
    axom::ArrayView<CoordType> m_x {}, m_y {};  //!< X,Y coordinates

    // New topology data.
    axom::ArrayView<ConnectivityType> m_connectivity {}, m_sizes {},
      m_offsets {};  //!< Connectivity data.

    // New field data.
    axom::ArrayView<axom::IndexType> m_original_zones {};  //!< View for originalZone field data.
    axom::ArrayView<double> m_norm_x {}, m_norm_y {};  //!< Fragment normals

    // New matset data
    axom::ArrayView<MaterialVF> m_volume_fractions {};
    axom::ArrayView<MaterialID> m_material_ids {}, m_mat_indices,
      m_mat_sizes {}, m_mat_offsets {};
  };

  /*!
   * \brief Return a view that can be copied to device.
   *
   * \return A view that be used to store data.
   */
  View view() { return m_view; }

private:
  View m_view;
};

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
    n_result["selectedZonesView"].set(selectedZonesView.data(),
                                      selectedZonesView.size());
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
                         const conduit::Node &AXOM_UNUSED_PARAM(n_fields),
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

    // Handle options.
    constexpr double DEFAULT_TOLERANCE = 1.e-10;
    constexpr int DEFAULT_MAX_ITERATIONS = 50;
    double tolerance = DEFAULT_TOLERANCE;
    if(n_options.has_child("tolerance") &&
       n_options["tolerance"].dtype().is_number())
    {
      tolerance = n_options["tolerance"].to_double();
    }
    int max_iterations = DEFAULT_MAX_ITERATIONS;
    if(n_options.has_child("max_iterations") &&
       n_options["max_iterations"].dtype().is_number())
    {
      max_iterations = n_options["max_iterations"].to_int();
    }

#ifdef AXOM_ELVIRA_DEBUG
    // Let's output the normals
    conduit::Node *n_result = new conduit::Node;
#endif

    // Count the number of fragments we'll make for the mixed zones.
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
#ifdef AXOM_ELVIRA_DEBUG
    conduit::Node &n_group1 = n_result->operator[]("group1");
    n_group1["mixedZones"].set(mixedZonesView.data(), mixedZonesView.size());
    n_group1["matZone"].set(matZoneView.data(), matZoneView.size());
    n_group1["matCount"].set(matCountView.data(), matCountView.size());
    n_group1["numFragments"].set(numFragments);
#endif
    AXOM_ANNOTATE_END("counting");

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
#ifdef AXOM_ELVIRA_DEBUG
    conduit::Node &n_group2 = n_result->operator[]("group2");
    n_group2["matZone"].set(matZoneView.data(), matZoneView.size());
    n_group2["matCount"].set(matCountView.data(), matCountView.size());
    n_group2["offsets"].set(matOffsetView.data(), matOffsetView.size());
#endif
    AXOM_ANNOTATE_END("offsets");

    //--------------------------------------------------------------------------
    // NOTE: This makes zone centers for all zones. That's okay since we use
    //       these values in making stencils.
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
#ifdef AXOM_ELVIRA_DEBUG
    conduit::Node &n_group3 = n_result->operator[]("group3");
    n_group3["xview"].set(xview.data(), xview.size());
    n_group3["yview"].set(yview.data(), yview.size());
    n_group3["zview"].set(zview.data(), zview.size());
#endif
    AXOM_ANNOTATE_END("centroids");

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

        // Reverse sort the materials by the volume fraction so the larger VFs are first.
        SLIC_ASSERT(ids.size() == matCount);
        axom::utilities::reverse_sort_multiple(vfs.data(), ids.data(), matCount);

        // Save sorted ids in sortedMaterialIdsView.
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          const auto fragmentIndex = offset + m;
          sortedMaterialIdsView[fragmentIndex] = ids[m];
        }

        // Retrieve the stencil data from neighbor zones.
        auto logical =
          deviceTopologyView.indexing().IndexToLogicalIndex(zoneIndex);
        for(int si = 0; si < StencilSize; si++)
        {
          // Stencil neighbor logical index.
          typename TopologyView::LogicalIndex neighbor(logical);

          // Neighbor offsets are in (-1, 0, 1) that get added to current zone's logical coordinate.
          const int neighborOffset[3] = {(si % 3) - 1,
                                         ((si % 9) / 3) - 1,
                                         (si / 9) - 1};
          for(int d = 0; d < ndims; d++)
          {
            neighbor[d] += neighborOffset[d];
          }

          // Clamp the neighbor to a zone that is inside the indexing space.
          neighbor = deviceTopologyView.indexing().clamp(neighbor);
          const auto neighborIndex = static_cast<typename MatsetView::ZoneIndex>(
            deviceTopologyView.indexing().LogicalIndexToIndex(neighbor));

          // Copy material vfs into the stencil.
          for(axom::IndexType m = 0; m < matCount; m++)
          {
            // Ask the neighbor zone for vf.
            const auto fragmentIndex = offset + m;
            typename MatsetView::FloatType vf = 0;
            matsetView.zoneContainsMaterial(neighborIndex,
                                            sortedMaterialIdsView[fragmentIndex],
                                            vf);

            // Store the vf into the stencil for the current material.
            const auto destIndex = fragmentIndex * StencilSize + si;
            fragmentVFStencilView[destIndex] = static_cast<double>(vf);
          }

          // The destination index for this coordinate stencil zone.
          const auto coordIndex = szIndex * StencilSize + si;

          // coord stencil
          xcStencilView[coordIndex] = xview[neighborIndex];
          ycStencilView[coordIndex] = yview[neighborIndex];
          zcStencilView[coordIndex] = zview.empty() ? 1. : zview[neighborIndex];
        }
      });
    AXOM_ANNOTATE_END("stencil");

    //--------------------------------------------------------------------------
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
        const double *fragmentVFStencilStart =
          fragmentVFStencilView.data() + offset * StencilSize;
        double *fragmentVectorsStart =
          fragmentVectorsView.data() + offset * numVectorComponents;

        // Produce normal for each material in this zone.
        int iskip = matCount - 1;
        elvira::elvira<ndims>::execute(matCount,
                                       fragmentVFStencilStart,
                                       fragmentVectorsStart,
                                       iskip);

#if defined(AXOM_ELVIRA_DEBUG) && !defined(AXOM_DEVICE_CODE)
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
#endif

        // Transform the normals.
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          double *normal =
            fragmentVectorsView.data() + ((offset + m) * numVectorComponents);
          elvira::transform(normal, jac);
#if defined(AXOM_ELVIRA_DEBUG) && !defined(AXOM_DEVICE_CODE)
          n_thisZone["mats"][m]["transformed_normal"].set(normal, 3);
#endif
        }
      });
#if defined(AXOM_ELVIRA_DEBUG)
    conduit::relay::io::save(*n_result, "elvira.yaml", "yaml");
    delete n_result;
#endif
    AXOM_ANNOTATE_END("vectors");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("allocate");

    // Free stencil data we no longer need.
    xcStencil.clear();
    ycStencil.clear();
    zcStencil.clear();

    // Determine the output type from the clip operations. Those are the shape
    // types that we're emitting into the MIR output. Create the builder.
    using CoordType = typename CoordsetView::value_type;
    using ClipResultType = typename std::conditional<
      ndims == 2,
      axom::primal::Polygon<CoordType, 2, axom::primal::PolygonArray::Static>,
      axom::primal::Polyhedron<CoordType, 3>>::type;
    using Builder =
      TopologyBuilder<ExecSpace, CoordsetView, TopologyView, MatsetView, ClipResultType, ndims>;
    Builder build;
    build.allocate(numFragments, n_newCoordset, n_newTopo, n_newFields, n_newMatset);
    if(n_matset.has_path("material_map"))
    {
      n_newMatset["material_map"].set(n_matset["material_map"]);
    }
    auto buildView = build.view();
    AXOM_ANNOTATE_END("allocate");
    //--------------------------------------------------------------------------

    AXOM_ANNOTATE_BEGIN("fragments");
    using VectorType = axom::primal::Vector<CoordType, ndims>;
    using PointType = axom::primal::Point<CoordType, ndims>;
    using PlaneType = axom::primal::Plane<CoordType, ndims>;
    using ShapeView = bputils::PrimalAdaptor<TopologyView, CoordsetView>;
    const ShapeView deviceShapeView {m_topologyView, m_coordsetView};
    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        const auto matCount = matCountView[szIndex];
        // Where this zone's fragment data starts.
        const auto offset = matOffsetView[szIndex];

        // Get the starting shape.
        auto shape = deviceShapeView.getShape(zoneIndex);

        // Get the zone's actual volume.
        double zoneVol = bputils::ComputeShapeAmount<ndims>::execute(shape);

        // NOTE: When we make fragments, we want to make the biggest ones first.
        //       With that in mind, we should probably reverse sort the materials based on VF above.
        //       If I do it there, I don't have to think about going backwards in loops, etc.
        for(axom::IndexType m = 0; m < matCount - 1; m++)
        {
          const auto fragmentIndex = offset + m;
          // Get this material fragment's normal and material id.
          const auto matId = sortedMaterialIdsView[fragmentIndex];
          const double *normalPtr =
            fragmentVectorsView.data() + (fragmentIndex * numVectorComponents);

          // Compute the desired fragment volume.
          // Get current material vf from the stencil. (should be faster than material view)
          constexpr int StencilCenter =
            (ndims == 3) ? 13 : ((ndims == 2) ? 4 : 1);
          const auto si = fragmentIndex * StencilSize + StencilCenter;
          const auto matVolume = zoneVol * fragmentVFStencilView[si];

          // Compute the shape bounding box.
          const auto bbox = axom::primal::compute_bounding_box(shape);

          // Make a plane situated at the bbox origin.
          VectorType normal(normalPtr, ndims);
          PlaneType P(normal, bbox.getCentroid(), false);
          //std::cout << "===================================================================\n";
          //std::cout << "Zone: " << zoneIndex << ", m=" << m << ", matId=" << matId << std::endl;
          //std::cout << "P=" << P << std::endl;
          // Compute distances from all points in shape to plane.
          PointType range[2];
          range[0] = range[1] = bbox.getCentroid();
          //std::cout << "range[0]=" << range[0] << std::endl;
          //std::cout << "range[1]=" << range[1] << std::endl;
          double dist[2] = {0., 0.};
          for(axom::IndexType ip = 0; ip < shape.numVertices(); ip++)
          {
            double d = P.signedDistance(shape[ip]);
            //std::cout << "ip=" << ip << ", shape[ip]=" << shape[ip] << ", " << d << std::endl;
            if(d < dist[0])
            {
              dist[0] = d;
              range[0] = bbox.getCentroid() + (normal * d);
              //std::cout << "dist[0]=" << d << ", range[0]=" << range[0] << std::endl;
            }
            if(d > dist[1])
            {
              dist[1] = d;
              range[1] = bbox.getCentroid() + (normal * d);
              //std::cout << "dist[1]=" << d << ", range[1]=" << range[1] << std::endl;
            }
          }

          bool searching = true;
          int iterations = 0;
          while(searching)
          {
            // Pick the middle of the range and position the plane there.
            const auto pt = PointType::lerp(range[0], range[1], 0.5);

            // The ELVIRA normals point away from the material. Axom's clipping
            // keeps the shape where the normal points into the shape. Reverse
            // the ELVIRA normal when forming the plane so we keep the right piece.
            P = PlaneType(-normal, pt, false);

            // Clip the shape at the current plane.
            auto clippedShape = axom::primal::clip(shape, P);

            // Find the volume of the clipped shape.
            const double fragmentVolume =
              bputils::ComputeShapeAmount<ndims>::execute(clippedShape);
            const double volumeError =
              axom::utilities::abs(matVolume - fragmentVolume);
            //std::cout << "\titerations=" << iterations << ", P=" << P << ", clippedShape=" << clippedShape << ", matVolume=" << matVolume << ", fragmentVolume=" << fragmentVolume << ", volumeError=" << volumeError << std::endl;
            if((volumeError <= tolerance) || (iterations >= max_iterations))
            {
              searching = false;

              // Emit clippedShape as material matId
              //std::cout << "zone=" << zoneIndex << ", fragmentIndex=" << fragmentIndex << ", mat=" << matId << ", iterations=" << iterations << ", clipped=" << clippedShape << std::endl;
              buildView.addShape(zoneIndex,
                                 fragmentIndex,
                                 clippedShape,
                                 matId,
                                 normalPtr);

              // Clip in the other direction to get the remaining fragment for the next material.
              P = PlaneType(normal, pt, false);
              shape = axom::primal::clip(shape, P);
            }
            else if(fragmentVolume < matVolume)
            {
              range[0] = pt;
            }
            else
            {
              range[1] = pt;
            }
            iterations++;
          }
        }

        // Emit the last leftover fragment.
        {
          const auto fragmentIndex = offset + matCount - 1;
          const auto matId = sortedMaterialIdsView[fragmentIndex];
          const double *normalPtr =
            fragmentVectorsView.data() + (fragmentIndex * numVectorComponents);
          //std::cout << "zone=" << zoneIndex << ", fragmentIndex=" << fragmentIndex << ", mat=" << matId << ", clipped=" << shape << std::endl;
          buildView.addShape(zoneIndex, fragmentIndex, shape, matId, normalPtr);
        }
      });
    AXOM_ANNOTATE_END("fragments");

    //--------------------------------------------------------------------------
#if defined(AXOM_ELVIRA_DEBUG)
    AXOM_ANNOTATE_BEGIN("save");
    conduit::Node n_mesh;
    n_mesh[n_newCoordset.path()].set_external(n_newCoordset);
    n_mesh[n_newTopo.path()].set_external(n_newTopo);
    n_mesh[n_newFields.path()].set_external(n_newFields);
    n_mesh[n_newMatset.path()].set_external(n_newMatset);

    // Verify the MIR output.
    conduit::Node info;
    if(!conduit::blueprint::mesh::verify(n_mesh, info))
    {
      info.print();
    }

    // Save the MIR output.
    conduit::relay::io::save(n_mesh, "elvira_output.yaml", "yaml");
    conduit::relay::io::blueprint::save_mesh(n_mesh, "elvira_output", "hdf5");
    AXOM_ANNOTATE_END("save");
#endif
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif
