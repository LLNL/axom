// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_TOPOLOGY_MAPPER_HPP_
#define AXOM_MIR_TOPOLOGY_MAPPER_HPP_

#include "axom/config.hpp"
#include "axom/CLI11.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/spin.hpp"

#include "axom/mir/utilities/PrimalAdaptor.hpp"
#include "axom/mir/utilities/VariableShape.hpp"
#include "axom/mir/utilities/utilities.hpp"

#include <conduit.hpp>
#include <conduit_blueprint.hpp>
#include <conduit_relay.hpp>

#include <iostream>

// Uncomment to emit debugging messages
//#define AXOM_DEBUG_TOPOLOGY_MAPPER

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
namespace detail
{

//------------------------------------------------------------------------------
/*!
 * \brief Return area where 2 polygons overlap.
 * \param shape1 The subject polygon.
 * \param shape2 The clip polygon.
 * \return The area common to 2 polygons.
 */
AXOM_SUPPRESS_HD_WARN
template <typename T, axom::primal::PolygonArray ARRAY_TYPE, int MAX_VERTS>
AXOM_HOST_DEVICE double shapeOverlap(
  const axom::primal::Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> &shape1,
  const axom::primal::Polygon<T, 2, ARRAY_TYPE, MAX_VERTS> &shape2,
  double eps = 1.e-10)
{
  constexpr bool tryFixOrientation = false;
  const auto p = axom::primal::clip(shape1, shape2, eps, tryFixOrientation);
  return p.area();
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Tetrahedron<T, 3> &shape1,
                                     const axom::primal::Tetrahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Hexahedron<T, 3> &shape1,
                                     const axom::primal::Hexahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Tetrahedron<T, 3> &shape1,
                                     const axom::primal::Hexahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const axom::primal::Hexahedron<T, 3> &shape1,
                                     const axom::primal::Tetrahedron<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const auto ph = axom::primal::clip(shape1, shape2, eps);
  return ph.volume();
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T, typename Shape2Type>
AXOM_HOST_DEVICE double shapeOverlap(const VariableShape<T, 3> &shape1,
                                     const Shape2Type &shape2,
                                     double eps = 1.e-10)
{
  const int id = shape1.id();
  double retval = 0.;
  if(id == axom::mir::views::Tet_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tet(shape1[0], shape1[1], shape1[2], shape1[3]);
    retval = shapeOverlap(tet, shape2, eps);
  }
  else if(id == axom::mir::views::Pyramid_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[2];
    shape1.splitPyramid(tets);
    for(int i = 0; i < 2; i++)
    {
      retval += shapeOverlap(tets[i], shape2, eps);
    }
  }
  else if(id == axom::mir::views::Wedge_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[3];
    shape1.splitWedge(tets);
    for(int i = 0; i < 3; i++)
    {
      const auto vol = shapeOverlap(tets[i], shape2, eps);
      retval += vol;
    }
  }
  else if(id == axom::mir::views::Hex_ShapeID)
  {
    axom::primal::Hexahedron<T, 3> hex(shape1[0],
                                       shape1[1],
                                       shape1[2],
                                       shape1[3],
                                       shape1[4],
                                       shape1[5],
                                       shape1[6],
                                       shape1[7]);
    retval = shapeOverlap(hex, shape2, eps);
  }
  else
  {
    assert("Unsupported shape type");
  }
  return retval;
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T, typename Shape1Type>
AXOM_HOST_DEVICE double shapeOverlap(const Shape1Type &shape1,
                                     const VariableShape<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  const int id = shape2.id();
  double retval = 0.;
  if(id == axom::mir::views::Tet_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tet(shape2[0], shape2[1], shape2[2], shape2[3]);
    retval = shapeOverlap(shape1, tet, eps);
  }
  else if(id == axom::mir::views::Pyramid_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[2];
    shape2.splitPyramid(tets);
    for(int i = 0; i < 2; i++)
    {
      retval += shapeOverlap(shape1, tets[i], eps);
    }
  }
  else if(id == axom::mir::views::Wedge_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[3];
    shape2.splitWedge(tets);
    for(int i = 0; i < 3; i++)
    {
      const auto vol = shapeOverlap(shape1, tets[i], eps);
      retval += vol;
    }
  }
  else if(id == axom::mir::views::Hex_ShapeID)
  {
    axom::primal::Hexahedron<T, 3> hex(shape2[0],
                                       shape2[1],
                                       shape2[2],
                                       shape2[3],
                                       shape2[4],
                                       shape2[5],
                                       shape2[6],
                                       shape2[7]);
    retval = shapeOverlap(shape1, hex, eps);
  }
  else
  {
    assert("Unsupported shape type");
  }
  return retval;
}

/*!
 * \brief Return the volume of the overlap between the shapes.
 * \param shape1 The subject shape.
 * \param shape2 The clip shape.
 * \return The volume of the overlap between the shapes.
 */
template <typename T>
AXOM_HOST_DEVICE double shapeOverlap(const VariableShape<T, 3> &shape1,
                                     const VariableShape<T, 3> &shape2,
                                     double eps = 1.e-10)
{
  int id = shape1.id();
  double retval = 0.;
  if(id == axom::mir::views::Tet_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tet(shape2[0], shape2[1], shape2[2], shape2[3]);
    retval = shapeOverlap(shape1, tet, eps);
  }
  else if(id == axom::mir::views::Pyramid_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[2];
    shape2.splitPyramid(tets);
    for(int i = 0; i < 2; i++)
    {
      retval += shapeOverlap(shape1, tets[i], eps);
    }
  }
  else if(id == axom::mir::views::Wedge_ShapeID)
  {
    axom::primal::Tetrahedron<T, 3> tets[3];
    shape2.splitWedge(tets);
    for(int i = 0; i < 3; i++)
    {
      const auto vol = shapeOverlap(shape1, tets[i], eps);
      retval += vol;
    }
  }
  else if(id == axom::mir::views::Hex_ShapeID)
  {
    axom::primal::Hexahedron<T, 3> hex(shape2[0],
                                       shape2[1],
                                       shape2[2],
                                       shape2[3],
                                       shape2[4],
                                       shape2[5],
                                       shape2[6],
                                       shape2[7]);
    retval = shapeOverlap(shape1, hex, eps);
  }
  else
  {
    assert("Unsupported shape type");
  }
  return retval;
}

}  // end namespace detail

//------------------------------------------------------------------------------
/**
 * \brief Take a source topology with clean matset and a target topology and
 *        intersect the source and target zones to build up a new matset on
 *        the target mesh that represents how the source mesh overlaps the
 *        target mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will execute.
 * \tparam SrcTopologyView The view type for the source topology.
 * \tparam SrcCoordsetView The view type for the source coordset.
 * \tparam TargetTopologyView The view type for the target topology.
 * \tparam TargetCoordsetView The view type for the target coordset.
 *
 * \note The use of topology and coordset views as template parameters allows
 *       this class to be instantiated for use with various topology and coordset
 *       types.
 */
template <typename ExecSpace,
          typename SrcTopologyView,
          typename SrcCoordsetView,
          typename TargetTopologyView,
          typename TargetCoordsetView>
class TopologyMapper
{
public:
  static_assert(SrcCoordsetView::dimension() == TargetCoordsetView::dimension(),
                "coordset dimension mismatch");

  using SrcShapeView = PrimalAdaptor<SrcTopologyView, SrcCoordsetView>;
  using TargetShapeView = PrimalAdaptor<TargetTopologyView, TargetCoordsetView>;

  /**
   * \brief Constructor
   *
   * \param srcTopoView The source topology view.
   * \param srcCoordsetView The source coordset view.
   * \param targetTopoView The target topology view.
   * \param targetCoordsetView The target coordset view.
   */
  TopologyMapper(const SrcTopologyView &srcTopoView,
                 const SrcCoordsetView &srcCoordsetView,
                 const TargetTopologyView &targetTopoView,
                 const TargetCoordsetView &targetCoordsetView)
    : m_srcView({srcTopoView, srcCoordsetView})
    , m_targetView({targetTopoView, targetCoordsetView})
  { }

  /**
   * \brief Intersect the source and target topologies and map the source
   *        material onto the target mesh.
   *
   * \param n_srcMesh The Conduit node that contains the coordset, topology, matset for the source mesh.
   * \param n_options A Conduit node that contains the algorithm options.
   * \param n_targetMesh The node that contains the coordset and topology for the source mesh.
   *
   * \verbatim
   * The n_options node must contain the following keys.
   *
   * source/matsetName: "name of source mesh matset"
   * target/topologyName: "name of target mesh topology"
   * target/matsetName "name of new target mesh matset"
   *
   * The n_options node may contain these optional parameters.
   *
   * source/selectedZones: [0,1,2,...]  # List of selected zone ids in source mesh.
   * target/selectedZones: [0,1,2,...]  # List of selected zone ids in target mesh.
   *
   * \endverbatim
   *
   * \note After executing, the n_targetMesh node will contain a new matset containing
   *       the results of the intersections with the src/target meshes.
   */
  void execute(const conduit::Node &n_srcMesh,
               const conduit::Node &n_options,
               conduit::Node &n_targetMesh)
  {
    AXOM_ANNOTATE_SCOPE("TopologyMapper::execute");

    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    const char *SRC_MATSET_NAME = "source/matsetName";
    const char *SRC_SELECTED_ZONES = "source/selectedZones";
    const char *TARGET_TOPOLOGY_NAME = "target/topologyName";
    const char *TARGET_MATSET_NAME = "target/matsetName";
    const char *TARGET_SELECTED_ZONES = "target/selectedZones";

    // Ensure required options exist.
    const char *required[] = {SRC_MATSET_NAME,
                              TARGET_TOPOLOGY_NAME,
                              TARGET_MATSET_NAME};
    for(const auto &key : required)
    {
      if(!n_options.has_path(key))
      {
        SLIC_ERROR(axom::fmt::format("Key \"{}\" missing from options.", key));
        return;
      }
    }
    const std::string srcMatsetName = n_options[SRC_MATSET_NAME].as_string();
    const std::string targetTopologyName =
      n_options[TARGET_TOPOLOGY_NAME].as_string();
    const std::string targetMatsetName =
      n_options[TARGET_MATSET_NAME].as_string();

    // Look at the source mesh's matset. Count the number of materials.
    const conduit::Node &n_matset =
      n_srcMesh.fetch_existing("matsets/" + srcMatsetName);
    const conduit::Node &n_materialMap =
      n_matset.fetch_existing("material_map");
    const auto nmats = n_materialMap.number_of_children();

    // Build up BVH that contains the src polygon bounding boxes.
    using SrcBoundingBox = typename SrcShapeView::BoundingBox;
    using src_value_type = typename SrcCoordsetView::value_type;
    AXOM_ANNOTATE_BEGIN("bbox");
    const auto srcView = m_srcView;
    SelectedZones<ExecSpace> srcSelection(srcView.numberOfZones(),
                                          n_options,
                                          SRC_SELECTED_ZONES);
    srcSelection.setSorted(false);
    const auto srcSelectionView = srcSelection.view();
    const axom::IndexType nSrcZones = srcSelectionView.size();
    axom::Array<SrcBoundingBox> srcBoundingBoxes(nSrcZones, nSrcZones, allocatorID);
    auto srcBoundingBoxesView = srcBoundingBoxes.view();
    axom::for_all<ExecSpace>(
      nSrcZones,
      AXOM_LAMBDA(axom::IndexType index) {
        const auto zi = srcSelectionView[index];
        srcBoundingBoxesView[index] = srcView.getBoundingBox(zi);
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
        std::cout << "source zone " << zi
                  << ": bbox=" << srcBoundingBoxesView[index] << std::endl;
#endif
      });
    AXOM_ANNOTATE_END("bbox");

    AXOM_ANNOTATE_BEGIN("build");
    axom::spin::BVH<SrcCoordsetView::dimension(), ExecSpace, src_value_type> bvh;
    bvh.setAllocatorID(allocatorID);
    bvh.initialize(srcBoundingBoxesView, srcBoundingBoxesView.size());
    AXOM_ANNOTATE_END("build");

    // -------------------------------------------------------------------------
    // Set up storage for a new matset.
    AXOM_ANNOTATE_BEGIN("allocation");
    ConduitAllocateThroughAxom<ExecSpace> c2a;

    conduit::Node &n_targetMatset = n_targetMesh["matsets/" + targetMatsetName];
    n_targetMatset["material_map"].set(n_materialMap);
    n_targetMatset["topology"].set(targetTopologyName);

    conduit::Node &n_volume_fractions = n_targetMatset["volume_fractions"];
    conduit::Node &n_material_ids = n_targetMatset["material_ids"];
    conduit::Node &n_indices = n_targetMatset["indices"];
    conduit::Node &n_sizes = n_targetMatset["sizes"];
    conduit::Node &n_offsets = n_targetMatset["offsets"];

    n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
    n_material_ids.set_allocator(c2a.getConduitAllocatorID());
    n_indices.set_allocator(c2a.getConduitAllocatorID());
    n_sizes.set_allocator(c2a.getConduitAllocatorID());
    n_offsets.set_allocator(c2a.getConduitAllocatorID());

    const axom::IndexType nTargetZones = m_targetView.numberOfZones();
    n_volume_fractions.set(conduit::DataType::float64(nmats * nTargetZones));
    n_material_ids.set(conduit::DataType::int32(nmats * nTargetZones));
    n_sizes.set(conduit::DataType::int32(nTargetZones));
    n_offsets.set(conduit::DataType::int32(nTargetZones));

    auto material_ids = make_array_view<int>(n_material_ids);
    auto volume_fractions = make_array_view<double>(n_volume_fractions);
    auto sizes = make_array_view<int>(n_sizes);
    auto offsets = make_array_view<int>(n_offsets);
    axom::mir::utilities::fill<ExecSpace>(volume_fractions, 0.);
    constexpr int MaterialEmpty = -1;
    axom::mir::utilities::fill<ExecSpace>(material_ids, MaterialEmpty);
    axom::mir::utilities::fill<ExecSpace>(sizes, 0);
    AXOM_ANNOTATE_END("allocation");

    // -------------------------------------------------------------------------
    // Iterate over the target zones and intersect them with source zones.
    AXOM_ANNOTATE_BEGIN("intersection");
    const auto targetView = m_targetView;
    RAJA::ReduceSum<reduce_policy, int> reduceSize(0);
    const auto srcMatIds =
      make_array_view<std::int64_t>(n_matset["material_ids"]);
    const auto bvh_device = bvh.getTraverser();
    SelectedZones<ExecSpace> targetSelection(targetView.numberOfZones(),
                                             n_options,
                                             TARGET_SELECTED_ZONES);
    targetSelection.setSorted(false);
    const auto targetSelectionView = targetSelection.view();
    axom::for_all<ExecSpace>(
      targetSelectionView.size(),
      AXOM_LAMBDA(axom::IndexType index) {
        // Get the target zone as a primal shape.
        const axom::IndexType zi = targetSelectionView[index];
        const auto targetBBox = targetView.getBoundingBox(zi);
        const auto targetShape = targetView.getShape(zi);
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
        std::cout << "-------------------------------\ntarget zone " << zi
                  << ": " << targetShape << ", bbox=" << targetBBox << std::endl;
#endif
        // Get the area or volume of the target shape (depends on the dimension).
        double targetAmount =
          ComputeShapeAmount<TargetCoordsetView::dimension()>::execute(
            targetShape);

        // Handle intersection in-depth of the bounding boxes intersected.
        auto handleIntersection = [&](std::int32_t currentNode,
                                      const std::int32_t *leafNodes) {
          const auto srcBboxIndex = leafNodes[currentNode];
          const auto srcZone = srcSelectionView[srcBboxIndex];
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
          std::cout << "handleIntersection: targetZone=" << zi
                    << ", srcZone=" << srcZone << std::endl;
#endif
          // Get the current zone as a primal shape.
          const auto srcShape = srcView.getShape(srcZone);

          // Determine polygon overlaps.
          constexpr double eps = 1.e-3;

          // Determine the overlap of the src and target shapes.
          double srcOverlapsTarget =
            detail::shapeOverlap(srcShape, targetShape, eps);

          if(srcOverlapsTarget > 0.)
          {
            double vf = srcOverlapsTarget / targetAmount;

            // Get the src material
            int mat = srcMatIds[srcZone];

#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
            std::cout << "\tintersection:" << std::endl
                      << "\t\ttargetShape=" << targetShape << std::endl
                      << "\t\tsrcShape=" << srcShape << std::endl
                      << "\t\tmat=" << mat << std::endl
                      << "\t\tsrcOverlapsTarget=" << srcOverlapsTarget
                      << std::endl
                      << "\t\ttargetAmount=" << targetAmount << std::endl
                      << "\t\tvf=" << vf << std::endl;
#endif

            // Add the src material into the target material.
            int *matids = material_ids.data() + zi * nmats;
            double *vfs = volume_fractions.data() + zi * nmats;
            for(int m = 0; m < nmats; m++)
            {
              if(matids[m] == mat)
              {
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
                std::cout << "\t\tAdded " << vf << " to slot " << m << std::endl;
#endif
                vfs[m] += vf;
                break;
              }
              else if(matids[m] == MaterialEmpty)
              {
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
                std::cout << "\t\tAdded new slot " << m << ", mat=" << mat
                          << ", vf=" << vf << std::endl;
#endif
                matids[m] = mat;
                vfs[m] = vf;
                sizes[zi]++;
                break;
              }
            }

            reduceSize += sizes[zi];
          }
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
          else
          {
            std::cout << "\tno intersection" << std::endl;
          }
#endif
        };

        // This predicate determines whether 2 bboxes intersect.
        auto bbIsect = [] /*AXOM_HOST_DEVICE*/ (
                         const SrcBoundingBox &queryBbox,
                         const SrcBoundingBox &bvhBbox) -> bool {
          bool rv = queryBbox.intersectsWith(bvhBbox);
#if defined(AXOM_DEBUG_TOPOLOGY_MAPPER) && !defined(AXOM_DEVICE_CODE)
          std::cout << "bbIsect: rv=" << rv << ", q=" << queryBbox
                    << ", bvh=" << bvhBbox << std::endl;
#endif
          return rv;
        };
        // Traverse BVH, looking for bboxes that intersect the current target bbox.
        bvh_device.traverse_tree(targetBBox, handleIntersection, bbIsect);
      });  // axom::for_all
    AXOM_ANNOTATE_END("intersection");

    // -------------------------------------------------------------------------
    // All the contributions have been added to the target matset. Finish building it.
    AXOM_ANNOTATE_BEGIN("finish");
    axom::exclusive_scan<ExecSpace>(sizes, offsets);
    const auto totalSize = reduceSize.get();
    n_indices.set(conduit::DataType::int32(totalSize));
    auto indices = make_array_view<int>(n_indices);
    axom::mir::utilities::fill<ExecSpace>(indices, -1);
    axom::for_all<ExecSpace>(
      targetView.numberOfZones(),
      AXOM_LAMBDA(axom::IndexType zi) {
        const auto start = offsets[zi];
        for(int i = 0; i < sizes[zi]; i++)
        {
          indices[start + i] = zi * nmats + i;
        }
      });
    AXOM_ANNOTATE_END("finish");
  }

  SrcShapeView m_srcView;
  TargetShapeView m_targetView;
};

}  // namespace blueprint
}  // namespace utilities
}  // namespace mir
}  // namespace axom

#endif
