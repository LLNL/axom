// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/primal/operators/intersection_volume.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/quest/detail/clipping/SphereClipper.hpp"
#include "axom/quest/detail/clipping/TetrahedronClipUtils.hpp"

namespace axom
{
namespace quest
{
namespace experimental
{

namespace
{

constexpr int SPHERE_MAX_SUBDIVISION_LEVELS = 4;
constexpr int SPHERE_BASE_REFINEMENT_LEVEL = 5;
constexpr double SPHERE_LINEARIZATION_EDGE_FRACTION = 0.08;

/*! \brief Precomputed sphere data used by the adaptive tet clipper. */
struct SphereClipData
{
  MeshClipperStrategy::Point3DType center;
  double radius;
  double radiusSquared;
  double targetEdgeSquared;
};

/*!
 * \brief Compute the local linearization scale for a sphere refinement level.
 *
 * \param [in] refinementLevel the sphere refinement level
 * \return the maximum edge length as a fraction of the sphere radius
 */
AXOM_MAYBE_UNUSED double sphereLinearizationEdgeFraction(int refinementLevel)
{
  double edgeFraction = SPHERE_LINEARIZATION_EDGE_FRACTION;
  for(int level = SPHERE_BASE_REFINEMENT_LEVEL; level < refinementLevel; ++level)
  {
    edgeFraction *= 0.5;
  }
  return edgeFraction;
}

/*!
 * \brief Choose the subdivision depth needed before local plane clipping.
 *
 * Coarse meshes need more refinement because the sphere curvature is visible
 * across a larger tet. Finer meshes can stop earlier and avoid unnecessary
 * subdivision work while still using the same GPU-friendly explicit stack.
 *
 * \param [in] tet the tetrahedron to examine
 * \param [in] sphere the sphere and target edge-length data
 * \return the number of eight-way subdivision levels
 */
AXOM_HOST_DEVICE
int chooseSubdivisionLevels(const MeshClipperStrategy::TetrahedronType& tet,
                            const SphereClipData& sphere)
{
  double edgeSq = detail::tetMaxEdgeSquared(tet);
  int levels = 0;
  while(levels < SPHERE_MAX_SUBDIVISION_LEVELS && edgeSq > sphere.targetEdgeSquared)
  {
    edgeSq *= 0.25;
    ++levels;
  }

  return levels;
}

/*!
 * \brief Approximate the sphere/tet overlap by linearizing the signed-distance
 * field over a small tetrahedron and clipping against the resulting plane.
 *
 * The signed-distance samples at the tet vertices define an affine field whose
 * zero isosurface is the best local planar approximation available from those
 * samples. The plane is oriented so that its positive side matches the sphere
 * interior, which matches primal::clip(tet, plane).
 *
 * \param [in] sphere the sphere data
 * \param [in] sqDistances squared distances from the sphere center to the tet vertices
 * \param [in] tetVolume the tetrahedron volume
 * \return the approximate sphere/tet overlap volume
 */
AXOM_HOST_DEVICE
double clipTetAgainstLinearizedSphere(const SphereClipData& sphere,
                                      const detail::TetVertexValues& sqDistances,
                                      double tetVolume)
{
  detail::TetVertexValues values;
  for(int i = 0; i < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++i)
  {
    values[i] = sphere.radius - sqrt(sqDistances[i]);
  }
  return detail::clipTetByVertexValues(values, tetVolume);
}

/*!
 * \brief Conservatively classify a tetrahedron against the sphere.
 *
 * This mirrors SphereClipper::polyhedronToLabel() but is specialized to tets
 * so the file-local adaptive clipper can use the same inexpensive screening.
 *
 * \param [in] tet the tetrahedron to classify
 * \param [in] sphere the sphere data
 * \param [out] sqDistances squared distances from the sphere center to the tet vertices
 * \return whether \p tet is inside, outside, or on the sphere boundary
 * \note The boundary classification is conservative.
 */
AXOM_HOST_DEVICE
MeshClipperStrategy::LabelType tetToSphereLabel(const MeshClipperStrategy::TetrahedronType& tet,
                                                const SphereClipData& sphere,
                                                detail::TetVertexValues& sqDistances)
{
  using LabelType = MeshClipperStrategy::LabelType;
  using Point3DType = MeshClipperStrategy::Point3DType;

  Point3DType minCoords = tet[0];
  Point3DType maxCoords = tet[0];
  bool allVerticesInside = true;

  for(int i = 0; i < MeshClipperStrategy::TetrahedronType::NUM_VERTS; ++i)
  {
    const auto& vert = tet[i];
    double sqDistToVert = 0.0;
    for(int dim = 0; dim < 3; ++dim)
    {
      minCoords[dim] = axom::utilities::min(minCoords[dim], vert[dim]);
      maxCoords[dim] = axom::utilities::max(maxCoords[dim], vert[dim]);

      const double delta = sphere.center[dim] - vert[dim];
      sqDistToVert += delta * delta;
    }
    sqDistances[i] = sqDistToVert;
    allVerticesInside = allVerticesInside && (sqDistToVert <= sphere.radiusSquared);
  }

  double sqDistToBb = 0.0;
  for(int dim = 0; dim < 3; ++dim)
  {
    const double delta = sphere.center[dim] < minCoords[dim]
      ? minCoords[dim] - sphere.center[dim]
      : (sphere.center[dim] > maxCoords[dim] ? sphere.center[dim] - maxCoords[dim] : 0.0);
    sqDistToBb += delta * delta;
  }

  if(sqDistToBb >= sphere.radiusSquared)
  {
    return LabelType::LABEL_OUT;
  }

  return allVerticesInside ? LabelType::LABEL_IN : LabelType::LABEL_ON;
}

/*!
 * \brief Compute a sphere/tet overlap using fixed-depth adaptive subdivision.
 *
 * Most boundary tets quickly break into fully interior or exterior subtets.
 * Only the smallest unresolved subtets pay for the local plane clip, which is
 * why this path is much cheaper than sending every boundary tet through the
 * generic sphere discretization and BVH clip pipeline.
 *
 * \param [in] tet the tetrahedron to clip
 * \param [in] sphere the sphere and target edge-length data
 * \param [in] rootVolume the volume of \p tet
 * \return the approximate sphere/tet overlap volume
 */
AXOM_MAYBE_UNUSED
AXOM_HOST_DEVICE
double clipTetAgainstSphere(const MeshClipperStrategy::TetrahedronType& tet,
                            const SphereClipData& sphere,
                            double rootVolume)
{
  using LabelType = MeshClipperStrategy::LabelType;
  using TetrahedronType = MeshClipperStrategy::TetrahedronType;
  constexpr int CHILD_COUNT = 8;
  constexpr int MAX_STACK_SIZE = 1 + (CHILD_COUNT - 1) * SPHERE_MAX_SUBDIVISION_LEVELS;

  const int maxSubdivisionLevels = chooseSubdivisionLevels(tet, sphere);
  double volume = 0.0;
  axom::StackArray<TetrahedronType, MAX_STACK_SIZE> tetStack;
  axom::StackArray<int, MAX_STACK_SIZE> levelStack;
  int stackSize = 1;
  tetStack[0] = tet;
  levelStack[0] = 0;
  axom::StackArray<double, SPHERE_MAX_SUBDIVISION_LEVELS + 1> levelVolumes;
  levelVolumes[0] = rootVolume;
  for(int level = 1; level <= maxSubdivisionLevels; ++level)
  {
    levelVolumes[level] = levelVolumes[level - 1] * 0.125;
  }

  // Use an explicit stack instead of recursion so the same code remains viable
  // in device execution spaces.
  while(stackSize > 0)
  {
    --stackSize;
    const TetrahedronType& currentTet = tetStack[stackSize];
    const int currentLevel = levelStack[stackSize];
    detail::TetVertexValues sqDistances;
    const LabelType label = tetToSphereLabel(currentTet, sphere, sqDistances);
    if(label == LabelType::LABEL_IN)
    {
      volume += levelVolumes[currentLevel];
      continue;
    }
    if(label == LabelType::LABEL_OUT)
    {
      continue;
    }

    if(currentLevel >= maxSubdivisionLevels)
    {
      volume += clipTetAgainstLinearizedSphere(sphere, sqDistances, levelVolumes[currentLevel]);
      continue;
    }

    detail::subdivideTetByMidpoints(currentTet, &tetStack[stackSize]);
    const int childLevel = currentLevel + 1;
    for(int i = 0; i < CHILD_COUNT; ++i)
    {
      // The stack bound is fixed by the subdivision depth and fanout above.
      levelStack[stackSize++] = childLevel;
    }
  }

  return volume;
}

}  // namespace

SphereClipper::SphereClipper(const klee::Geometry& kGeom, const std::string& name)
  : MeshClipperStrategy(kGeom)
  , m_name(name.empty() ? std::string("Sphere") : name)
  , m_transformer(m_extTrans)
{
  extractClipperInfo();

  transformSphere();
}

bool SphereClipper::labelCellsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                    axom::Array<LabelType>& labels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  int allocId = shapeMesh.getAllocatorID();
  auto cellCount = shapeMesh.getCellCount();
  if(labels.size() < cellCount || labels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    labels = axom::Array<LabelType>(ArrayOptions::Uninitialized(), cellCount, cellCount, allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelCellsInOutImpl<axom::SEQ_EXEC>(shapeMesh, labels.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelCellsInOutImpl<axom::OMP_EXEC>(shapeMesh, labels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelCellsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, labels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelCellsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, labels.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void SphereClipper::labelCellsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                        axom::ArrayView<LabelType> labels)
{
  auto cellCount = shapeMesh.getCellCount();
  auto cellsAsHexes = shapeMesh.getCellsAsHexes();
  auto cellVolumes = shapeMesh.getCellVolumes();
  constexpr double EPS = 1e-10;
  auto sphere = m_sphere;

  // Host policies use precomputed cell bounds and direct sphere tests.
  if constexpr(!axom::execution_space<ExecSpace>::onDevice())
  {
    AXOM_ANNOTATE_SCOPE("SphereClipper::label_cells_cpu");
    auto cellBounds = shapeMesh.getCellBoundingBoxes();
    const auto sphereCenter = sphere.getCenter();
    const double sphereRadiusSquared = sphere.getRadius() * sphere.getRadius();
    axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
      if(axom::utilities::isNearlyEqual(cellVolumes[cellId], 0.0, EPS) ||
         primal::squared_distance(sphereCenter, cellBounds[cellId]) >= sphereRadiusSquared)
      {
        labels[cellId] = LabelType::LABEL_OUT;
        return;
      }

      const auto& hex = cellsAsHexes[cellId];
      for(int vi = 0; vi < HexahedronType::NUM_HEX_VERTS; ++vi)
      {
        if(primal::squared_distance(sphereCenter, hex[vi]) > sphereRadiusSquared)
        {
          labels[cellId] = LabelType::LABEL_ON;
          return;
        }
      }
      labels[cellId] = LabelType::LABEL_IN;
    });
    return;
  }

  // Device policies use the general polyhedron classifier.
  AXOM_ANNOTATE_SCOPE("SphereClipper::label_cells_device");
  axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType cellId) {
    LabelType& cellLabel = labels[cellId];
    if(axom::utilities::isNearlyEqual(cellVolumes[cellId], 0.0, EPS))
    {
      cellLabel = LabelType::LABEL_OUT;
      return;
    }
    const auto& hex = cellsAsHexes[cellId];
    cellLabel = polyhedronToLabel(hex, sphere);
    // Note: cellLabel may be set to LABEL_ON if polyhedronToLabel
    // cannot efficiently determine whether the hex is IN or OUT.
    // See MeshClipperStrategy::labelCellsInOut().
  });
  return;
}

bool SphereClipper::labelTetsInOut(quest::experimental::ShapeMesh& shapeMesh,
                                   axom::ArrayView<const axom::IndexType> cellIds,
                                   axom::Array<LabelType>& tetLabels)
{
  SLIC_ERROR_IF(shapeMesh.dimension() != 3, "SphereClipper requires a 3D mesh.");

  const axom::IndexType cellCount = cellIds.size();
  const int allocId = shapeMesh.getAllocatorID();

  if(tetLabels.size() < cellCount * NUM_TETS_PER_HEX ||
     tetLabels.getAllocatorID() != shapeMesh.getAllocatorID())
  {
    tetLabels = axom::Array<LabelType>(ArrayOptions::Uninitialized(),
                                       cellCount * NUM_TETS_PER_HEX,
                                       cellCount * NUM_TETS_PER_HEX,
                                       allocId);
  }

  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    labelTetsInOutImpl<axom::SEQ_EXEC>(shapeMesh, cellIds, tetLabels.view());
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    labelTetsInOutImpl<axom::OMP_EXEC>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    labelTetsInOutImpl<axom::CUDA_EXEC<256>>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    labelTetsInOutImpl<axom::HIP_EXEC<256>>(shapeMesh, cellIds, tetLabels.view());
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return true;
}

template <typename ExecSpace>
void SphereClipper::labelTetsInOutImpl(quest::experimental::ShapeMesh& shapeMesh,
                                       axom::ArrayView<const axom::IndexType> cellIds,
                                       axom::ArrayView<LabelType> tetLabels)
{
  const axom::IndexType cellCount = cellIds.size();
  auto tetVolumes = shapeMesh.getTetVolumes();
  constexpr double EPS = 1e-10;
  auto sphere = m_sphere;

  // Host policies classify the mesh's precomputed tetrahedra directly.
  if constexpr(!axom::execution_space<ExecSpace>::onDevice())
  {
    AXOM_ANNOTATE_SCOPE("SphereClipper::label_tets_cpu");
    auto meshTets = shapeMesh.getCellsAsTets();
    axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType ci) {
      const axom::IndexType cellId = cellIds[ci];
      for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
      {
        LabelType& tetLabel = tetLabels[ci * NUM_TETS_PER_HEX + ti];
        const axom::IndexType tetId = cellId * NUM_TETS_PER_HEX + ti;
        if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
        {
          tetLabel = LabelType::LABEL_OUT;
          continue;
        }

        tetLabel = polyhedronToLabel(meshTets[tetId], sphere);
      }
    });
    return;
  }

  // Device policies use the hex decomposition and general tet classifier.
  AXOM_ANNOTATE_SCOPE("SphereClipper::label_tets_device");
  auto meshHexes = shapeMesh.getCellsAsHexes();

  axom::for_all<ExecSpace>(cellCount, [=] AXOM_HOST_DEVICE(axom::IndexType ci) {
    axom::IndexType cellId = cellIds[ci];
    const HexahedronType& hex = meshHexes[cellId];

    TetrahedronType cellTets[NUM_TETS_PER_HEX];
    ShapeMesh::hexToTets(hex, cellTets);

    for(IndexType ti = 0; ti < NUM_TETS_PER_HEX; ++ti)
    {
      LabelType& tetLabel = tetLabels[ci * NUM_TETS_PER_HEX + ti];
      const axom::IndexType tetId = cellId * NUM_TETS_PER_HEX + ti;
      if(axom::utilities::isNearlyEqual(tetVolumes[tetId], 0.0, EPS))
      {
        tetLabel = LabelType::LABEL_OUT;
        continue;
      }
      const TetrahedronType& tet = cellTets[ti];
      tetLabel = polyhedronToLabel(tet, sphere);
      // Note: cellLabel may be set to LABEL_ON if polyhedronToLabel
      // cannot efficiently determine whether the tet is IN or OUT.
      // See MeshClipperStrategy::labelTetsInOut().
    }
  });
  return;
}

template <typename Polyhedron>
AXOM_HOST_DEVICE inline MeshClipperStrategy::LabelType SphereClipper::polyhedronToLabel(
  const Polyhedron& verts,
  const SphereType& sphere)
{
  /*
    If bounding box of polyhedron is more than the radius distance
    from center, it is LABEL_OUT.  (Comparing vertices for this check
    can miss intersections by edges and facets, so we compare bounding
    box.)

    Otherwise, polyhedron is labeled either LABEL_ON or LABEL_IN.
    Sphere is convex, so polyhedron is IN only if all vertices are inside.

    Some polyhedra may be labeled ON even though they are actually outside.
    This is a conservative error that is fixed ty the clip function later.
    The purpose of labeling is bypass the clip function where we can do
    it efficiently.
  */
  BoundingBox3DType bb(verts[0]);
  auto vertCount = Polyhedron::numVertices();
  for(int i = 1; i < vertCount; ++i)
  {
    bb.addPoint(verts[i]);
  }

  const double sqRad = sphere.getRadius() * sphere.getRadius();

  double sqDistToBb = primal::squared_distance(sphere.getCenter(), bb);

  if(sqDistToBb >= sqRad)
  {
    return LabelType::LABEL_OUT;
  }

  for(int i = 0; i < vertCount; ++i)
  {
    const auto& vert = verts[i];
    double sqDistToVert = axom::primal::squared_distance(sphere.getCenter(), vert);
    if(sqDistToVert > sqRad)
    {
      return LabelType::LABEL_ON;
    }
  }
  return LabelType::LABEL_IN;
}

bool SphereClipper::getGeometryAsOcts(quest::experimental::ShapeMesh& shapeMesh,
                                      axom::Array<axom::primal::Octahedron<double, 3>>& octs)
{
  AXOM_ANNOTATE_SCOPE("SphereClipper::getGeometryAsOcts");
  int octCount = 0;
  axom::quest::discretize(m_sphereBeforeTrans, m_levelOfRefinement, octs, octCount);

  auto octsView = octs.view();
  auto transformer = m_transformer;
  const int allocId = shapeMesh.getAllocatorID();
  axom::for_all<axom::SEQ_EXEC>(octCount, [=] AXOM_HOST_DEVICE(axom::IndexType iOct) {
    OctahedronType& oct = octsView[iOct];
    for(int iVert = 0; iVert < OctType::NUM_VERTS; ++iVert)
    {
      Point3DType& ptCoords = oct[iVert];
      transformer.transform(ptCoords.array());
    }
  });

  // The discretize method uses host data. Place into the required allocator if needed.
  if(octs.getAllocatorID() != allocId)
  {
    octs = axom::Array<axom::primal::Octahedron<double, 3>>(octs, allocId);
  }

  SLIC_DEBUG(axom::fmt::format("SphereClipper '{}' {}-level refined got {} geometry octs.",
                               name(),
                               m_levelOfRefinement,
                               octs.size()));
  return true;
}

bool SphereClipper::specializedClipTets(quest::experimental::ShapeMesh& shapeMesh,
                                        axom::ArrayView<double> ovlap,
                                        const axom::ArrayView<IndexType>& tetIds,
                                        conduit::Node& statistics)
{
  // Dispatch the analytic sphere-tetrahedron overlap in the selected execution
  // space.
  bool didClip = false;
  switch(shapeMesh.getRuntimePolicy())
  {
  case axom::runtime_policy::Policy::seq:
    didClip = specializedClipTetsImpl<axom::SEQ_EXEC>(shapeMesh, ovlap, tetIds, statistics);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case axom::runtime_policy::Policy::omp:
    didClip = specializedClipTetsImpl<axom::OMP_EXEC>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case axom::runtime_policy::Policy::cuda:
    didClip = specializedClipTetsImpl<axom::CUDA_EXEC<256>>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case axom::runtime_policy::Policy::hip:
    didClip = specializedClipTetsImpl<axom::HIP_EXEC<256>>(shapeMesh, ovlap, tetIds, statistics);
    break;
#endif
  default:
    SLIC_ERROR("Axom Internal error: Unhandled execution policy.");
  }
  return didClip;
}

template <typename ExecSpace>
bool SphereClipper::specializedClipTetsImpl(quest::experimental::ShapeMesh& shapeMesh,
                                            axom::ArrayView<double> ovlap,
                                            const axom::ArrayView<IndexType>& tetIds,
                                            conduit::Node& statistics)
{
  AXOM_ANNOTATE_SCOPE("SphereClipper::specialized_clip");
  auto meshTets = shapeMesh.getCellsAsTets();
  const auto sphere = m_sphere;
  const IndexType tetCount = tetIds.size();

  AXOM_ANNOTATE_BEGIN("SphereClipper::exact_overlap");
  axom::for_all<ExecSpace>(tetCount, [=] AXOM_HOST_DEVICE(axom::IndexType ti) {
    const axom::IndexType tetId = tetIds[ti];
    const axom::IndexType cellId = tetId / NUM_TETS_PER_HEX;
    SLIC_ASSERT(cellId >= 0 && cellId < ovlap.size());
    const auto& tet = meshTets[tetId];
    const double volume = primal::intersection_volume(sphere, tet);
    detail::addToOverlapVolume<ExecSpace>(ovlap.data() + cellId, volume);
  });
  AXOM_ANNOTATE_END("SphereClipper::exact_overlap");

  statistics["onSum"].set_int64(tetCount);
  statistics["clipsSum"].set_int64(tetCount);
  return true;
}

void SphereClipper::extractClipperInfo()
{
  const auto c = m_info.fetch_existing("center").as_double_array();
  const double radius = m_info.fetch_existing("radius").as_double();
  Point3DType center;
  for(int d = 0; d < 3; ++d)
  {
    center[d] = c[d];
  }
  m_sphereBeforeTrans = SphereType(center, radius);
  m_levelOfRefinement = m_info.fetch_existing("levelOfRefinement").to_int32();
}

// Include external transformations in m_sphere.
void SphereClipper::transformSphere()
{
  const auto& centerBeforeTrans = m_sphereBeforeTrans.getCenter();
  const double radiusBeforeTrans = m_sphereBeforeTrans.getRadius();
  Point3DType surfacePtBeforeTrans {centerBeforeTrans.array() +
                                    Point3DType::NumericArray {radiusBeforeTrans, 0, 0}};

  auto center = m_transformer.getTransformed(centerBeforeTrans);
  Point3DType surfacePoint = m_transformer.getTransformed(surfacePtBeforeTrans);
  const double radius = Vector3DType(center, surfacePoint).norm();
  m_sphere = SphereType(center, radius);
}

}  // namespace experimental
}  // end namespace quest
}  // end namespace axom
