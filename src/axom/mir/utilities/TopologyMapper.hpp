// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_TOPOLOGY_MAPPER_HPP_
#define AXOM_MIR_TOPOLOGY_MAPPER_HPP_

#include <axom/config.hpp>
#include <axom/CLI11.hpp>
#include <axom/core.hpp>
#include <axom/mir.hpp>
#include <axom/primal.hpp>
#include <axom/slic.hpp>
#include <axom/spin.hpp>

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
/**
 * \brief Fill an ArrayView with a value.
 *
 * \tparam ExecSpace The execution space where the fill will be done.
 * \tparam T The data type of the values in the ArrayView.
 *
 * \param view The ArrayView being filled.
 * \param fillValue The value to be used for filling the ArrayView.
 */
template <typename ExecSpace, typename T>
void fill(axom::ArrayView<T> view, T fillValue)
{
  axom::for_all<ExecSpace>(
    view.size(),
    AXOM_LAMBDA(axom::IndexType index) { view[index] = fillValue; });
}

//------------------------------------------------------------------------------
/**
 * \brief Represent various shapes that consist of points.
 *
 * \tparam T The coordinate precision type.
 * \tparam NDIMS The spatial coordinate dimension.
 * \tparam N The max number of points that the shape can contain.
 */
template <typename T, int NDIMS, int N = 8>
class VariableShape
{
public:
  using PointType = axom::primal::Point<T, NDIMS>;

  /*!
   * \brief Return the shape id type.
   * \return The shape id type.
   */
  AXOM_HOST_DEVICE int id() const { return m_shapeId; }

  /*!
   * \brief Return the number of points in the shape.
   * \return The number of points in the shape.
   */
  AXOM_HOST_DEVICE axom::IndexType size() const { return m_points.size(); }

  /*!
   * \brief Add a point to the shape.
   * \param pt The point to add.
   */
  AXOM_HOST_DEVICE void push_back(const PointType &pt)
  {
    m_points.push_back(pt);
  }

  /*!
   * \brief Return the \a index'th point.
   * \param index The index of the point to return.
   * \return The desired point.
   */
  AXOM_HOST_DEVICE const PointType &operator[](axom::IndexType index) const
  {
    return m_points[index];
  }

  /*!
   * \brief Return unsigned shape volume.
   * \return Unsigned shape volume.
   */
  AXOM_HOST_DEVICE double volume() const
  {
    double retval = 0.;
    if(m_shapeId == axom::mir::views::Tet_ShapeID)
    {
      axom::primal::Tetrahedron<T, 3> tet(m_points[0],
                                          m_points[1],
                                          m_points[2],
                                          m_points[3]);
      retval = tet.volume();
    }
    else if(m_shapeId == axom::mir::views::Pyramid_ShapeID)
    {
      axom::primal::Tetrahedron<T, 3> tets[2];
      splitPyramid(tets);
      for(int i = 0; i < 2; i++)
      {
        retval += tets[i].volume();
      }
    }
    else if(m_shapeId == axom::mir::views::Wedge_ShapeID)
    {
      axom::primal::Tetrahedron<T, 3> tets[3];
      splitWedge(tets);
      for(int i = 0; i < 3; i++)
      {
        retval += tets[i].volume();
      }
    }
    else if(m_shapeId == axom::mir::views::Hex_ShapeID)
    {
      axom::primal::Hexahedron<T, 3> hex(m_points[0],
                                         m_points[1],
                                         m_points[2],
                                         m_points[3],
                                         m_points[4],
                                         m_points[5],
                                         m_points[6],
                                         m_points[7]);
      retval = hex.volume();
    }
    else
    {
      assert("Unsupported shape type");
    }
    return retval;
  }

  /*!
   * \brief Split the shape into tets.
   * \param[out] tets The output array of tets that make up the shape.
   */
  AXOM_HOST_DEVICE void splitPyramid(axom::primal::Tetrahedron<T, NDIMS> tets[2]) const
  {
    assert(m_shapeId == axom::mir::views::Pyramid_ShapeID);
    tets[0] = axom::primal::Tetrahedron<T, NDIMS>(m_points[0],
                                                  m_points[1],
                                                  m_points[3],
                                                  m_points[4]);
    tets[1] = axom::primal::Tetrahedron<T, NDIMS>(m_points[1],
                                                  m_points[2],
                                                  m_points[3],
                                                  m_points[4]);
  }

  /*!
   * \brief Split the shape into tets.
   * \param[out] tets The output array of tets that make up the shape.
   */
  AXOM_HOST_DEVICE void splitWedge(axom::primal::Tetrahedron<T, NDIMS> tets[3]) const
  {
    assert(m_shapeId == axom::mir::views::Wedge_ShapeID);
    tets[0] = axom::primal::Tetrahedron<T, NDIMS>(m_points[0],
                                                  m_points[1],
                                                  m_points[2],
                                                  m_points[3]);
    tets[1] = axom::primal::Tetrahedron<T, NDIMS>(m_points[3],
                                                  m_points[1],
                                                  m_points[2],
                                                  m_points[5]);
    tets[2] = axom::primal::Tetrahedron<T, NDIMS>(m_points[3],
                                                  m_points[4],
                                                  m_points[1],
                                                  m_points[5]);
  }

  int m_shapeId;
  axom::StaticArray<PointType, N> m_points;
};

template <typename T, int NDIMS, int N = 8>
std::ostream &operator<<(std::ostream &os, const VariableShape<T, NDIMS, N> &obj)
{
  os << "{shapeId=" << obj.m_shapeId << ", points={";
  for(int i = 0; i < obj.m_points.size(); i++)
  {
    if(i > 0)
    {
      os << ", ";
    }
    os << obj.m_points[i];
  }
  os << "}}";
  return os;
}

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

#if defined(DIFFERENT_CLIPPER)
  using Polygon = axom::primal::Polygon<T, 2, ARRAY_TYPE, MAX_VERTS>;
  Polygon overlapPoly;
  Clipper<Polygon>::Clip(shape1, shape2, overlapPoly, eps);
  return overlapPoly.area();
#else
  const auto p = axom::primal::clip(shape1, shape2, eps, tryFixOrientation);
  return p.area();
#endif
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

/*!
 * \brief Base template for computing a shape's area or volume.
 */
template <int NDIMS>
struct ComputeShapeAmount
{ };

/*!
 * \brief 2D specialization for shapes to compute area.
 */
template <>
struct ComputeShapeAmount<2>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.area();
  }
};

/*!
 * \brief 3D specialization for shapes to compute volume.
 */
template <>
struct ComputeShapeAmount<3>
{
  template <typename ShapeType>
  static inline AXOM_HOST_DEVICE double execute(const ShapeType &shape)
  {
    return shape.volume();
  }
};

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

  /**
   * \brief This class is a view that exposes shapes in MIR topology and
   *        coordset views as primal shapes. SFINAE is used for the zone()
   *        method so various shapes can be returned.
   */
  template <typename TopologyView, typename CoordsetView>
  struct PrimalAdaptor
  {
    using value_type = typename CoordsetView::value_type;
    using Polygon = axom::primal::Polygon<value_type,
                                          CoordsetView::dimension(),
                                          axom::primal::PolygonArray::Static>;
    using Tetrahedron =
      axom::primal::Tetrahedron<value_type, CoordsetView::dimension()>;
    using Hexahedron =
      axom::primal::Hexahedron<value_type, CoordsetView::dimension()>;
    using BoundingBox =
      axom::primal::BoundingBox<value_type, CoordsetView::dimension()>;

    /**
     * \brief Return the number of zones in the associated topology view.
     * \return The number of zones in the associated topology view.
     */
    AXOM_HOST_DEVICE axom::IndexType numZones() const
    {
      return m_topoView.numberOfZones();
    }

    /**
     * \brief Return the bounding box for the zi'th zone.
     *
     * \param zi The index of the zone whose bounding box we want.
     *
     * \return The bounding box of the zi'th zone.
     */
    AXOM_HOST_DEVICE BoundingBox getBoundingBox(axom::IndexType zi) const
    {
      const auto zone = m_topoView.zone(zi);
      const auto nnodes = zone.numberOfNodes();
      BoundingBox b;
      for(axom::IndexType i = 0; i < nnodes; i++)
      {
        b.addPoint(m_coordsetView[zone.getId(i)]);
      }
      return b;
    }

    // SFINAE methods for returning a zone as a primal shape (or VariableShape
    // if primal lacks the shape type).

    /**
     * \brief Get the zone \a zi as a polygon. This is enabled when the input topology is 2D.
     *
     * \param zi The index of the zone we want to return.
     *
     * \return A Polygon that represents the zone from the input topology.
     */
    template <int TDIM = CoordsetView::dimension()>
    AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, Polygon>::type getShape(
      axom::IndexType zi) const
    {
      const auto zone = m_topoView.zone(zi);
      Polygon p;
      for(axom::IndexType i = 0; i < zone.numberOfNodes(); i++)
      {
        p.addVertex(m_coordsetView[zone.getId(i)]);
      }
      return p;
    }

    /**
     * \brief Get the zone \a zi as a Tetrahedron. This is enabled when the input topology contains tets.
     *
     * \param zi The index of the zone we want to return.
     *
     * \return A Tetrahedron that represents the zone from the input topology.
     */
    template <int TDIM = CoordsetView::dimension(),
              typename ShapeType = typename TopologyView::ShapeType>
    AXOM_HOST_DEVICE typename std::enable_if<
      TDIM == 3 &&
        std::is_same<ShapeType,
                     axom::mir::views::TetShape<typename ShapeType::ConnectivityStorage>>::value,
      Tetrahedron>::type
    getShape(axom::IndexType zi) const
    {
      const auto zone = m_topoView.zone(zi);
      return Tetrahedron(m_coordsetView[zone.getId(0)],
                         m_coordsetView[zone.getId(1)],
                         m_coordsetView[zone.getId(2)],
                         m_coordsetView[zone.getId(3)]);
    }

    /**
     * \brief Get the zone \a zi as a Hexahedron. This is enabled when the input topology contains hexs.
     *
     * \param zi The index of the zone we want to return.
     *
     * \return A Hexahedron that represents the zone from the input topology.
     */
    template <int TDIM = CoordsetView::dimension(),
              typename ShapeType = typename TopologyView::ShapeType>
    AXOM_HOST_DEVICE typename std::enable_if<
      TDIM == 3 &&
        std::is_same<ShapeType,
                     axom::mir::views::HexShape<typename ShapeType::ConnectivityStorage>>::value,
      Hexahedron>::type
    getShape(axom::IndexType zi) const
    {
      const auto zone = m_topoView.zone(zi);
      return Hexahedron(m_coordsetView[zone.getId(0)],
                        m_coordsetView[zone.getId(1)],
                        m_coordsetView[zone.getId(2)],
                        m_coordsetView[zone.getId(3)],
                        m_coordsetView[zone.getId(4)],
                        m_coordsetView[zone.getId(5)],
                        m_coordsetView[zone.getId(6)],
                        m_coordsetView[zone.getId(7)]);
    }

    /**
     * \brief Get the zone \a zi as a VariableShape. This is enabled when the input
     *        topology contains wedges, pyramids, or mixed shapes.
     *
     * \param zi The index of the zone we want to return.
     *
     * \return A VariableShape that represents the zone from the input topology.
     *         VariableShape is used because primal lacks Pyramid, Wedge classes and
     *         because if we're using this, we might have a mesh with mixed zone types.
     */
    template <int TDIM = CoordsetView::dimension(),
              typename ShapeType = typename TopologyView::ShapeType>
    AXOM_HOST_DEVICE typename std::enable_if<
      (TDIM == 3) &&
        (std::is_same<ShapeType,
                      axom::mir::views::PyramidShape<typename ShapeType::ConnectivityType>>::value ||
         std::is_same<ShapeType,
                      axom::mir::views::WedgeShape<typename ShapeType::ConnectivityType>>::value ||
         std::is_same<ShapeType,
                      axom::mir::views::VariableShape<
                        typename ShapeType::ConnectivityType>>::value),
      detail::VariableShape<value_type, 3>>::type
    getShape(axom::IndexType zi) const
    {
      const auto zone = m_topoView.zone(zi);
      detail::VariableShape<value_type, 3> shape;
      shape.m_shapeId = zone.id();
      for(int i = 0; i < zone.numberOfNodes(); i++)
      {
        shape.push_back(m_coordsetView[zone.getId(i)]);
      }
      return shape;
    }

    TopologyView m_topoView;
    CoordsetView m_coordsetView;
  };

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
    namespace bputils = axom::mir::utilities::blueprint;
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
    bputils::SelectedZones<ExecSpace> srcSelection(srcView.numZones(),
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
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

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

    const axom::IndexType nTargetZones = m_targetView.numZones();
    n_volume_fractions.set(conduit::DataType::float64(nmats * nTargetZones));
    n_material_ids.set(conduit::DataType::int32(nmats * nTargetZones));
    n_sizes.set(conduit::DataType::int32(nTargetZones));
    n_offsets.set(conduit::DataType::int32(nTargetZones));

    auto material_ids = bputils::make_array_view<int>(n_material_ids);
    auto volume_fractions = bputils::make_array_view<double>(n_volume_fractions);
    auto sizes = bputils::make_array_view<int>(n_sizes);
    auto offsets = bputils::make_array_view<int>(n_offsets);
    detail::fill<ExecSpace>(volume_fractions, 0.);
    constexpr int MaterialEmpty = -1;
    detail::fill<ExecSpace>(material_ids, MaterialEmpty);
    detail::fill<ExecSpace>(sizes, 0);
    AXOM_ANNOTATE_END("allocation");

    // -------------------------------------------------------------------------
    // Iterate over the target zones and intersect them with source zones.
    AXOM_ANNOTATE_BEGIN("intersection");
    const auto targetView = m_targetView;
    RAJA::ReduceSum<reduce_policy, int> reduceSize(0);
    const auto srcMatIds =
      bputils::make_array_view<std::int64_t>(n_matset["material_ids"]);
    const auto bvh_device = bvh.getTraverser();
    bputils::SelectedZones<ExecSpace> targetSelection(targetView.numZones(),
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
          detail::ComputeShapeAmount<TargetCoordsetView::dimension()>::execute(
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
        auto bbIsect = [] AXOM_HOST_DEVICE(const SrcBoundingBox &queryBbox,
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
    auto indices = bputils::make_array_view<int>(n_indices);
    detail::fill<ExecSpace>(indices, -1);
    axom::for_all<ExecSpace>(
      targetView.numZones(),
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
