// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file IntersectionShaper.hpp
 *
 * \brief Helper class for intersection-based shaping queries
 */

#ifndef AXOM_QUEST_INTERSECTION_SHAPER__HPP_
#define AXOM_QUEST_INTERSECTION_SHAPER__HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"
#include "axom/klee.hpp"

#ifndef AXOM_USE_MFEM
  #error Shaping functionality requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif

#include "axom/quest/Shaper.hpp"
#include "axom/spin/BVH.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/detail/shaping/shaping_helpers.hpp"

#include "mfem.hpp"

#include "axom/fmt.hpp"
#include "axom/fmt/locale.h"

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

// clang-format off
#if defined (AXOM_USE_RAJA) && defined (AXOM_USE_UMPIRE)
  using seq_exec = axom::SEQ_EXEC;

  #if defined(AXOM_USE_OPENMP)
    using omp_exec = axom::OMP_EXEC;
  #else
    using omp_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_CUDA)
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
  #else
    using cuda_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_HIP)
    constexpr int HIP_BLOCK_SIZE = 256;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
  #else
    using hip_exec = seq_exec;
  #endif
#endif
// clang-format on

namespace axom
{
namespace quest
{
class IntersectionShaper : public Shaper
{
public:
  using BoundingBoxType = primal::BoundingBox<double, 3>;
  using PolyhedronType = primal::Polyhedron<double, 3>;
  using OctahedronType = primal::Octahedron<double, 3>;
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using TetrahedronType = primal::Tetrahedron<double, 3>;

  /// Choose runtime policy for RAJA
  enum ExecPolicy
  {
    seq = 0,
    omp = 1,
    cuda = 2,
    hip = 3
  };

public:
  IntersectionShaper(const klee::ShapeSet& shapeSet,
                     sidre::MFEMSidreDataCollection* dc)
    : Shaper(shapeSet, dc)
  { }

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  virtual ~IntersectionShaper()
  {
    deallocateReplacementRuleStorage();
  }
#endif

  //@{
  //!  @name Functions to get and set shaping parameters related to intersection; supplements parameters in base class

  void setLevel(int level) { m_level = level; }

  void setExecPolicy(int policy) { m_execPolicy = (ExecPolicy)policy; }
  //@}

private:
  /**
   * \brief Helper method to decompose a hexahedron Polyhedron
   *        into 24 Tetrahedrons.
   *        Each Tetrahedron consists of 4 points:
   *          (1) The mean of all Polyhedron points (centroid)
   *          (2,3) Two adjacent vertices on a Polyhedron face
   *          (4) The mean of the current Polyhedron face
   *
   * \param poly [in] The hexahedron Polyhedron to decompose
   * \param tets [out] The tetrahedrons
   *
   * \note Assumes given Polyhedron has 8 points
   * \note Assumes tets is pre-allocated
   */

  AXOM_HOST_DEVICE
  void decompose_hex_to_tets(const PolyhedronType& poly, TetrahedronType* tets)
  {
    // Hex center (hc)
    Point3D hc = poly.centroid();

    //Face means (fm)
    Point3D fm1 = Point3D::midpoint(Point3D::midpoint(poly[0], poly[1]),
                                    Point3D::midpoint(poly[2], poly[3]));

    Point3D fm2 = Point3D::midpoint(Point3D::midpoint(poly[0], poly[1]),
                                    Point3D::midpoint(poly[4], poly[5]));

    Point3D fm3 = Point3D::midpoint(Point3D::midpoint(poly[0], poly[3]),
                                    Point3D::midpoint(poly[4], poly[7]));

    Point3D fm4 = Point3D::midpoint(Point3D::midpoint(poly[1], poly[2]),
                                    Point3D::midpoint(poly[5], poly[6]));

    Point3D fm5 = Point3D::midpoint(Point3D::midpoint(poly[2], poly[3]),
                                    Point3D::midpoint(poly[6], poly[7]));

    Point3D fm6 = Point3D::midpoint(Point3D::midpoint(poly[4], poly[5]),
                                    Point3D::midpoint(poly[6], poly[7]));

    // Initialize tets
    tets[0] = TetrahedronType(hc, poly[1], poly[0], fm1);
    tets[1] = TetrahedronType(hc, poly[0], poly[3], fm1);
    tets[2] = TetrahedronType(hc, poly[3], poly[2], fm1);
    tets[3] = TetrahedronType(hc, poly[2], poly[1], fm1);

    tets[4] = TetrahedronType(hc, poly[4], poly[0], fm2);
    tets[5] = TetrahedronType(hc, poly[0], poly[1], fm2);
    tets[6] = TetrahedronType(hc, poly[1], poly[5], fm2);
    tets[7] = TetrahedronType(hc, poly[5], poly[4], fm2);

    tets[8] = TetrahedronType(hc, poly[3], poly[0], fm3);
    tets[9] = TetrahedronType(hc, poly[0], poly[4], fm3);
    tets[10] = TetrahedronType(hc, poly[4], poly[7], fm3);
    tets[11] = TetrahedronType(hc, poly[7], poly[3], fm3);

    tets[12] = TetrahedronType(hc, poly[5], poly[1], fm4);
    tets[13] = TetrahedronType(hc, poly[1], poly[2], fm4);
    tets[14] = TetrahedronType(hc, poly[2], poly[6], fm4);
    tets[15] = TetrahedronType(hc, poly[6], poly[5], fm4);

    tets[16] = TetrahedronType(hc, poly[6], poly[2], fm5);
    tets[17] = TetrahedronType(hc, poly[2], poly[3], fm5);
    tets[18] = TetrahedronType(hc, poly[3], poly[7], fm5);
    tets[19] = TetrahedronType(hc, poly[7], poly[6], fm5);

    tets[20] = TetrahedronType(hc, poly[7], poly[4], fm6);
    tets[21] = TetrahedronType(hc, poly[4], poly[5], fm6);
    tets[22] = TetrahedronType(hc, poly[5], poly[6], fm6);
    tets[23] = TetrahedronType(hc, poly[6], poly[7], fm6);
  }

  /**
   * \brief Helper method to check if an Octahedron has duplicate
   *        vertices
   *
   * \param poly [in] The Octahedron to check for duplicate vertices
   * \return True if duplicate vertices found, false otherwise
   */
  AXOM_HOST_DEVICE
  bool oct_has_duplicate_verts(const OctahedronType& oct) const
  {
    for(int i = 0; i < 6; i++)
    {
      for(int j = i + 1; j < OctahedronType::NUM_OCT_VERTS; j++)
      {
        // operator= for Point does not want to play nice...
        if(oct[i][0] == oct[j][0] && oct[i][1] == oct[j][1] &&
           oct[i][2] == oct[j][2])
        {
          return true;
        }
      }
    }
    return false;
  }

public:
//@{
//!  @name Functions related to the stages for a given shape

/// Initializes the spatial index for shaping
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  template <typename ExecSpace>
  void prepareShapeQueryImpl(klee::Dimensions shapeDimension,
                             const klee::Shape& shape)
  {
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(
        "Running intersection-based shaper in execution Space: {}",
        axom::execution_space<ExecSpace>::name())));

    // Save current/default allocator
    const int current_allocator = axom::getDefaultAllocatorID();

    // Determine new allocator (for CUDA/HIP policy, set to Unified)
    // Set new default to device
    axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

    const auto& shapeName = shape.getName();
    AXOM_UNUSED_VAR(shapeDimension);
    AXOM_UNUSED_VAR(shapeName);

    SLIC_INFO(axom::fmt::format("Current shape is {}", shapeName));

    // Number of points in polyline
    int pointcount = getSurfaceMesh()->getNumberOfNodes();

    Point2D* polyline = axom::allocate<Point2D>(pointcount);

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(" Refinement level set to {} ", m_level)));

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(
        " Checking contour with {} points for degenerate segments ",
        pointcount)));

    enum
    {
      R = 1,
      Z = 0
    };

    // Add contour points
    int polyline_size = 0;
    const double EPS = m_vertexWeldThreshold;
    const double EPS_SQ = EPS * EPS;
    for(int i = 0; i < pointcount; ++i)
    {
      Point2D cur_point;
      m_surfaceMesh->getNode(i, cur_point.data());

      // Check for degenerate segments
      if(polyline_size > 0)
      {
        using axom::utilities::isNearlyEqual;
        const Point2D& prev_point = polyline[polyline_size - 1];

        // both r are (nearly) equal to 0
        if(isNearlyEqual(cur_point[R], 0.0, EPS) &&
           isNearlyEqual(prev_point[R], 0.0, EPS))
        {
          continue;
        }

        // both Z are (nearly) the same
        if(isNearlyEqual(cur_point[Z], prev_point[Z], EPS))
        {
          continue;
        }

        // both points are (nearly) the same
        if(squared_distance(cur_point, prev_point) < EPS_SQ)
        {
          continue;
        }
      }

      polyline[polyline_size] = cur_point;
      polyline_size += 1;
    }

    // Check if we need to flip the points order.
    // discretize() is only valid if x increases as index increases
    bool flip = false;
    if(polyline_size > 1)
    {
      if(polyline[0][Z] > polyline[1][Z])
      {
        flip = true;
        SLIC_INFO("Order of contour points has been reversed!"
                  << " Discretization algorithm expects Z values of contour "
                     "points to be increasing");
      }
    }

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(" Discretizing contour with {} points ", polyline_size)));

    // Flip point order
    if(flip)
    {
      int i = polyline_size - 1;
      int j = 0;
      while(i > j)
      {
        axom::utilities::swap<Point2D>(polyline[i], polyline[j]);
        i -= 1;
        j += 1;
      }
    }

    // Generate the Octahedra
    const bool disc_status = axom::quest::discretize<ExecSpace>(polyline,
                                                                polyline_size,
                                                                m_level,
                                                                m_octs,
                                                                m_octcount);

    // Oddities required by hip to avoid capturing `this`
    OctahedronType* local_octs = m_octs;

    AXOM_UNUSED_VAR(disc_status);  // silence warnings in release configs
    SLIC_ASSERT_MSG(
      disc_status,
      "Discretization of contour has failed. Check that contour is valid");

    SLIC_INFO(
      axom::fmt::format("Contour has been discretized into {} octahedra ",
                        m_octcount));

    if(this->isVerbose())
    {
      // Print out the bounding box containing all the octahedra
      BoundingBoxType all_oct_bb;
      for(int i = 0; i < m_octcount; i++)
      {
        all_oct_bb.addBox(primal::compute_bounding_box(m_octs[i]));
      }
      SLIC_INFO(axom::fmt::format(
        "DEBUG: Bounding box containing all generated octahedra "
        "has dimensions:\n\t{}",
        all_oct_bb));

      // Print out the total volume of all the octahedra
      using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
      RAJA::ReduceSum<REDUCE_POL, double> total_oct_vol(0.0);
      axom::for_all<ExecSpace>(
        m_octcount,
        AXOM_LAMBDA(axom::IndexType i) {
          // Convert Octahedron into Polyhedrom
          PolyhedronType octPoly;
          double octVolume;

          octPoly.addVertex(local_octs[i][0]);
          octPoly.addVertex(local_octs[i][1]);
          octPoly.addVertex(local_octs[i][2]);
          octPoly.addVertex(local_octs[i][3]);
          octPoly.addVertex(local_octs[i][4]);
          octPoly.addVertex(local_octs[i][5]);

          octPoly.addNeighbors(0, {1, 5, 4, 2});
          octPoly.addNeighbors(1, {0, 2, 3, 5});
          octPoly.addNeighbors(2, {0, 4, 3, 1});
          octPoly.addNeighbors(3, {1, 2, 4, 5});
          octPoly.addNeighbors(4, {0, 5, 3, 2});
          octPoly.addNeighbors(5, {0, 1, 3, 4});

          octVolume = octPoly.volume();

          // Flip sign if volume is negative
          // (expected when vertex order is reversed)
          if(octVolume < 0)
          {
            octVolume = -octVolume;
          }
          total_oct_vol += octVolume;
        });

      SLIC_INFO(axom::fmt::format(
        "DEBUG: Total volume of all generated octahedra is {}",
        total_oct_vol.get()));

      // Check if any Octahedron are degenerate with all points {0,0,0}
      RAJA::ReduceSum<REDUCE_POL, int> num_degenerate(0);
      axom::for_all<ExecSpace>(
        m_octcount,
        AXOM_LAMBDA(axom::IndexType i) {
          OctahedronType degenerate_oct;
          if(local_octs[i].equals(degenerate_oct))
          {
            num_degenerate += 1;
          }
        });

      SLIC_INFO(
        axom::fmt::format("DEBUG: {} Octahedron found with all points (0,0,0)",
                          num_degenerate.get()));

      // Dump discretized octs as a tet mesh
      axom::mint::Mesh* tetmesh;
      axom::quest::mesh_from_discretized_polyline(m_octs,
                                                  m_octcount,
                                                  polyline_size - 1,
                                                  tetmesh);
      axom::mint::write_vtk(tetmesh, "discretized_surface_of_revolution.vtk");
      delete tetmesh;

    }  // end of verbose output for contour

    axom::deallocate(polyline);
    axom::setDefaultAllocator(current_allocator);
  }
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  template <typename ExecSpace>
  void runShapeQueryImpl(const klee::Shape& shape)
  {
    constexpr int NUM_VERTS_PER_HEX = 8;
    constexpr int NUM_COMPS_PER_VERT = 3;
    constexpr int NUM_TETS_PER_HEX = 24;
    constexpr double ZERO_THRESHOLD = 1.e-10;

    // Save current/default allocator
    const int current_allocator = axom::getDefaultAllocatorID();

    // Determine new allocator (for CUDA/HIP policy, set to Unified)
    // Set new default to device
    axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

    int* ZERO = axom::allocate<int>(
      1,
      axom::getUmpireResourceAllocatorID(umpire::resource::Host));
    ZERO[0] = 0;

    SLIC_INFO(
      axom::fmt::format("{:-^80}",
                        " Inserting Octahedra bounding boxes into BVH "));

    // Generate the BVH tree over the octahedra
    // Access-aligned bounding boxes
    m_aabbs = axom::allocate<BoundingBoxType>(m_octcount);

    // Oddities required by hip to avoid capturing `this`
    OctahedronType* local_octs = m_octs;
    BoundingBoxType* local_aabbs = m_aabbs;

    // Get the bounding boxes for the Octahedrons
    axom::for_all<ExecSpace>(
      m_octcount,
      AXOM_LAMBDA(axom::IndexType i) {
        local_aabbs[i] = primal::compute_bounding_box<double, 3>(local_octs[i]);
      });

    // Insert Octahedra Bounding Boxes into BVH.
    //bvh.setAllocatorID(poolID);
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(m_aabbs, m_octcount);

    SLIC_INFO(axom::fmt::format("{:-^80}", " Querying the BVH tree "));

    mfem::Mesh* mesh = getDC()->GetMesh();

    // Intersection algorithm only works on linear elements
    SLIC_ASSERT(mesh != nullptr);
    int const NE = mesh->GetNE();
    m_num_elements = NE;

    if(this->isVerbose())
    {
      SLIC_INFO(axom::fmt::format(
        "{:-^80}",
        axom::fmt::format(
          " Initializing {} hexahedral elements from given mesh ",
          m_num_elements)));
    }

    if(NE > 0)
    {
      SLIC_ASSERT(mesh->GetNodes() == nullptr ||
                  mesh->GetNodes()->FESpace()->GetOrder(0));
    }

    // Create and register a scalar field for this shape's volume fractions
    // The Degrees of Freedom will be in correspondence with the elements
    auto* volFrac = this->newVolFracGridFunction();
    auto volFracName = axom::fmt::format("shape_vol_frac_{}", shape.getName());
    this->getDC()->RegisterField(volFracName, volFrac);

    // Initialize hexahedral elements
    m_hexes = axom::allocate<PolyhedronType>(NE);
    m_hex_bbs = axom::allocate<BoundingBoxType>(NE);

    // Oddities required by hip to avoid capturing `this`
    PolyhedronType* local_hexes = m_hexes;
    BoundingBoxType* local_hex_bbs = m_hex_bbs;

    // Initialize vertices from mfem mesh and
    // set each shape volume fraction to 1
    // Allocation size is:
    // # of elements * # of vertices per hex * # of components per vertex
    double* vertCoords =
      axom::allocate<double>(NE * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT);
    for(int i = 0; i < NE; i++)
    {
      // Get the indices of this element's vertices
      mfem::Array<int> verts;
      mesh->GetElementVertices(i, verts);
      SLIC_ASSERT(verts.Size() == NUM_VERTS_PER_HEX);

      // Set each shape volume fraction to 1
      double vf = 1.0;
      (*volFrac)(i) = vf;

      // Get the coordinates for the vertices
      for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
      {
        for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
        {
          vertCoords[(i * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT) +
                     (j * NUM_COMPS_PER_VERT) + k] =
            (mesh->GetVertex(verts[j]))[k];
        }
      }
    }

    // Initialize each hexahedral element and its bounding box
    axom::for_all<ExecSpace>(
      NE,
      AXOM_LAMBDA(axom::IndexType i) {
        // Set each hexahedral element vertices
        local_hexes[i] = PolyhedronType();
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIndex = (i * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT) +
            j * NUM_COMPS_PER_VERT;
          local_hexes[i].addVertex({vertCoords[vertIndex],
                                    vertCoords[vertIndex + 1],
                                    vertCoords[vertIndex + 2]});

          // Set hexahedra components to zero if within threshold
          if(axom::utilities::isNearlyEqual(local_hexes[i][j][0],
                                            0.0,
                                            ZERO_THRESHOLD))
          {
            local_hexes[i][j][0] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_hexes[i][j][1],
                                            0.0,
                                            ZERO_THRESHOLD))
          {
            local_hexes[i][j][1] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_hexes[i][j][2],
                                            0.0,
                                            ZERO_THRESHOLD))
          {
            local_hexes[i][j][2] = 0.0;
          }
        }

        // Get bounding box for hexahedral element
        local_hex_bbs[i] =
          primal::compute_bounding_box<double, 3>(local_hexes[i]);
      });  // end of loop to initialize hexahedral elements and bounding boxes

    // Deallocate no longer needed
    axom::deallocate(vertCoords);

    // Set octahedra components to zero if within threshold
    axom::for_all<ExecSpace>(
      m_octcount,
      AXOM_LAMBDA(axom::IndexType i) {
        for(int j = 0; j < OctahedronType::NUM_OCT_VERTS; j++)
        {
          if(axom::utilities::isNearlyEqual(local_octs[i][j][0],
                                            0.0,
                                            ZERO_THRESHOLD))
          {
            local_octs[i][j][0] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_octs[i][j][1],
                                            0.0,
                                            ZERO_THRESHOLD))
          {
            local_octs[i][j][1] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_octs[i][j][2],
                                            0.0,
                                            ZERO_THRESHOLD))
          {
            local_octs[i][j][2] = 0.0;
          }
        }
      });

    // Find which octahedra bounding boxes intersect hexahedron bounding boxes
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Finding octahedra candidates for each hexahedral element "));

    axom::Array<IndexType> offsets(NE);
    axom::Array<IndexType> counts(NE);
    axom::Array<IndexType> candidates;
    bvh.findBoundingBoxes(offsets,
                          counts,
                          candidates,
                          NE,
                          reinterpret_cast<BoundingBoxType*>(m_hex_bbs));

    //Deallocate no longer needed variables
    axom::deallocate(m_aabbs);

    // Get the total number of candidates
    using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
    using ATOMIC_POL = typename axom::execution_space<ExecSpace>::atomic_policy;

    const auto counts_v = counts.view();
    RAJA::ReduceSum<REDUCE_POL, int> totalCandidates(0);
    axom::for_all<ExecSpace>(
      NE,
      AXOM_LAMBDA(axom::IndexType i) { totalCandidates += counts_v[i]; });

    // Initialize hexahedron indices and octahedra candidates
    axom::IndexType* hexIndices =
      axom::allocate<axom::IndexType>(totalCandidates.get() * NUM_TETS_PER_HEX);
    axom::IndexType* octCandidates =
      axom::allocate<axom::IndexType>(totalCandidates.get() * NUM_TETS_PER_HEX);

    // Tetrahedrons from hexes (24 for each hex)
    TetrahedronType* tets =
      axom::allocate<TetrahedronType>(NE * NUM_TETS_PER_HEX);

    // Index into 'tets'
    axom::IndexType* tetIndices =
      axom::allocate<axom::IndexType>(totalCandidates.get() * NUM_TETS_PER_HEX);

    // New total number of candidates after omitting degenerate octahedra
    int* newTotalCandidates = axom::allocate<int>(1);
    axom::copy(newTotalCandidates, ZERO, sizeof(int));

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Decomposing each hexahedron element into 24 tetrahedrons "));

    AXOM_PERF_MARK_SECTION("init_tets",
                           axom::for_all<ExecSpace>(
                             NE,
                             AXOM_LAMBDA(axom::IndexType i) {
                               TetrahedronType cur_tets[NUM_TETS_PER_HEX];
                               decompose_hex_to_tets(local_hexes[i], cur_tets);

                               for(int j = 0; j < NUM_TETS_PER_HEX; j++)
                               {
                                 tets[i * NUM_TETS_PER_HEX + j] = cur_tets[j];
                               }
                             }););

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Linearizing each tetrahedron, octahedron candidate pair "));

    const auto offsets_v = offsets.view();
    const auto candidates_v = candidates.view();
    AXOM_PERF_MARK_SECTION(
      "init_candidates",
      axom::for_all<ExecSpace>(
        NE,
        AXOM_LAMBDA(axom::IndexType i) {
          for(int j = 0; j < counts_v[i]; j++)
          {
            int octIdx = candidates_v[offsets_v[i] + j];
            if(!oct_has_duplicate_verts(local_octs[octIdx]))
            {
              for(int k = 0; k < NUM_TETS_PER_HEX; k++)
              {
                auto idx = RAJA::atomicAdd<ATOMIC_POL>(newTotalCandidates, 1);
                hexIndices[idx] = i;
                octCandidates[idx] = octIdx;
                tetIndices[idx] = i * NUM_TETS_PER_HEX + k;
              }
            }
          }
        }););

    // Overlap volume is the volume of clip(oct,tet)
    m_overlap_volumes = axom::allocate<double>(NE);

    // Hex volume is the volume of the hexahedron element
    m_hex_volumes = axom::allocate<double>(NE);

    // Oddities required by hip to avoid capturing `this`
    double* local_overlap_volumes = m_overlap_volumes;
    double* local_hex_volumes = m_hex_volumes;

    // Set initial values to 0
    axom::for_all<ExecSpace>(
      NE,
      AXOM_LAMBDA(axom::IndexType i) {
        local_hex_volumes[i] = 0;
        local_overlap_volumes[i] = 0;
      });

    SLIC_INFO(
      axom::fmt::format("{:-^80}", " Calculating hexahedron element volume "));

    AXOM_PERF_MARK_SECTION("hex_volume",
                           axom::for_all<ExecSpace>(
                             NE * NUM_TETS_PER_HEX,
                             AXOM_LAMBDA(axom::IndexType i) {
                               double tet_volume = tets[i].volume();
                               RAJA::atomicAdd<ATOMIC_POL>(
                                 local_hex_volumes + (i / NUM_TETS_PER_HEX),
                                 tet_volume);
                             }););

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      " Calculating element overlap volume from each tet-oct pair "));

    AXOM_PERF_MARK_SECTION(
      "oct_tet_volume",
      axom::for_all<ExecSpace>(
        newTotalCandidates[0],
        AXOM_LAMBDA(axom::IndexType i) {
          int index = hexIndices[i];
          int octIndex = octCandidates[i];
          int tetIndex = tetIndices[i];
          PolyhedronType poly =
            primal::clip(local_octs[octIndex], tets[tetIndex]);

          // Poly is valid
          if(poly.numVertices() >= 4)
          {
            double clip_volume = poly.volume();
            // Flip sign if negative
            if(clip_volume < 0)
            {
              clip_volume = -clip_volume;
            }
            RAJA::atomicAdd<ATOMIC_POL>(local_overlap_volumes + index,
                                        clip_volume);
          }
        }););

    RAJA::ReduceSum<REDUCE_POL, double> totalOverlap(0);
    RAJA::ReduceSum<REDUCE_POL, double> totalHex(0);

    axom::for_all<ExecSpace>(
      NE,
      AXOM_LAMBDA(axom::IndexType i) {
        totalOverlap += local_overlap_volumes[i];
        totalHex += local_hex_volumes[i];
      });

    SLIC_INFO(axom::fmt::format("Total overlap volume with shape is {}",
                                this->allReduceSum(totalOverlap)));
    SLIC_INFO(axom::fmt::format("Total mesh volume is {}",
                                this->allReduceSum(totalHex)));

    // Deallocate no longer needed variables
    axom::deallocate(ZERO);
    axom::deallocate(m_hexes);
    axom::deallocate(m_hex_bbs);
    axom::deallocate(hexIndices);
    axom::deallocate(octCandidates);
    axom::deallocate(tets);
    axom::deallocate(tetIndices);
    axom::deallocate(newTotalCandidates);
    axom::deallocate(m_octs);

    axom::setDefaultAllocator(current_allocator);
  }  // end of runShapeQuery() function
#endif

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  // These methods are private in support of replacement rules.
private:

  std::string materialNameToFieldName(const std::string &materialName) const
  {
    const std::string vol_frac_fmt("vol_frac_{}");
    auto name = axom::fmt::format(vol_frac_fmt, materialName);
    return name;
  }

  std::string fieldNameToMaterialName(const std::string &fieldName) const
  {
    const std::string vol_frac_("vol_frac_");
    std::string name;
    if(fieldName.find(vol_frac_) == 0)
      name = fieldName.substr(vol_frac_.size());
    return name;
  }

  // Gets the grid function and material number for a material name.
  std::pair<mfem::GridFunction *, int> getMaterial(const std::string &materialName)
  {
    // If we already know about the material, return it.
    for(size_t i = 0; i < m_vf_material_names.size(); i++)
    {
      if(m_vf_material_names[i] == materialName)
      {
        return std::make_pair(m_vf_grid_functions[i], static_cast<int>(i));
      }
    }

    // Get or create the volume fraction field for this shape's material
    auto materialVolFracName = materialNameToFieldName(materialName);
    mfem::GridFunction* matVolFrac = nullptr;
    if(this->getDC()->HasField(materialVolFracName))
    {
      matVolFrac = this->getDC()->GetField(materialVolFracName);
std::cout << "*** Found existing GF for " << materialName << std::endl;
    }
    else
    {
std::cout << "*** Made new GF for " << materialName << std::endl;
      matVolFrac = newVolFracGridFunction();
      this->getDC()->RegisterField(materialVolFracName, matVolFrac);
      // Zero out the volume fractions (on host).
      memset(matVolFrac->begin(), 0, matVolFrac->Size() * sizeof(double));
    }

    // Add the material to our vectors.
    int idx = static_cast<int>(m_vf_grid_functions.size());
    m_vf_grid_functions.push_back(matVolFrac);
    m_vf_material_names.push_back(materialName);

    return std::make_pair(matVolFrac, idx);
  }

  // If we are passed in a data collection that already has volume fractions
  // in it, then we need to add them to our vectors. We maintain out own
  // vectors because we assume that the order of materials does not change.
  // The mfem::DataCollection uses a map internally so if we add materials,
  // we could change the traversal order.
  void populateMaterials()
  {
    std::vector<std::string> materialNames;
    for(auto it : this->getDC()->GetFieldMap())
    {
      std::string materialName = fieldNameToMaterialName(it.first);
      if(!materialName.empty())
        materialNames.emplace_back(materialName);
    }
    // Add any of these existing fields to this class' bookkeeping.
    for(const auto &materialName : materialNames)
    {
      (void)getMaterial(materialName);
    }
#if 1
    std::cout << "*** populateMaterials: names={";
    for(const auto &name : m_vf_material_names)
      std::cout << name << ", ";
    std::cout << "}" << std::endl;
#endif
  }

  template <typename ExecSpace>
  bool allocateReplacementRuleStorage(int size)
  {
    // Allocate the memory if we have not already.
    bool first = m_vf_sums == nullptr;
    if(first)
    {
      // Save current/default allocator
      const int current_allocator = axom::getDefaultAllocatorID();

      // Determine new allocator (for CUDA/HIP policy, set to Unified)
      // Set new default to device
      m_replacement_allocator = axom::execution_space<ExecSpace>::allocatorID();
      axom::setDefaultAllocator(m_replacement_allocator);

      m_vf_sums = axom::allocate<double>(size);
      m_vf_subtract = axom::allocate<double>(size);
      m_vf_occupied = axom::allocate<double>(size);

      axom::setDefaultAllocator(current_allocator);
    }

    return first;
  }

  void deallocateReplacementRuleStorage()
  {
    if(m_vf_sums)
    {
      const int current_allocator = axom::getDefaultAllocatorID();
      axom::setDefaultAllocator(m_replacement_allocator);
      axom::deallocate(m_vf_sums);
      m_vf_sums = nullptr;
      axom::deallocate(m_vf_subtract);
      m_vf_subtract = nullptr;
      axom::deallocate(m_vf_occupied);
      m_vf_occupied = nullptr;
      axom::setDefaultAllocator(current_allocator);
    }
  }

  // Ares-style replacement rules built using forall.
  template <typename ExecSpace>
  void applyReplacementRulesImpl(const klee::Shape& shape)
  {
    // Make sure the material lists are up to date.
    populateMaterials();

    // Get this shape's material, creating it if needed.
    auto matVF = getMaterial(shape.getMaterial());
    int dataSize = matVF.first->Size();

    // Allocate some memory for the replacement rule data arrays.
    bool first = allocateReplacementRuleStorage<ExecSpace>(dataSize);
int values_per_line = 20;
    // Determine which grid functions need to be considered for VF updates.
    std::vector<std::pair<mfem::GridFunction *,int>> gf_order_by_matnumber;
    std::vector<mfem::GridFunction *> updateVFs, excludeVFs;
    if(!shape.getMaterialsReplaced().empty())
    {
      // Include materials replaced in updateVFs.
      for(const auto &name : shape.getMaterialsReplaced())
      {
        gf_order_by_matnumber.emplace_back(getMaterial(name));
      }
    }
    else
    {
std::cout << "*** excludeVFs" << std::endl;
      // Include all materials except those in "does_not_replace".
      // We'll also sort them by material number since the field map
      // sorts them by name rather than order added.
      for(auto it : this->getDC()->GetFieldMap())
      {
        // Check whether the field name looks like a VF field.
        std::string name = fieldNameToMaterialName(it.first);
        if(!name.empty())
        {
          // See if the field is in the exclusion list. For the normal
          // case, the list is empty so we'd add the material.
          auto it2 = std::find(shape.getMaterialsNotReplaced().cbegin(),
                               shape.getMaterialsNotReplaced().cend(),
                               name);
          // The field is not in the exclusion list so add it to vfs.
          if(it2 == shape.getMaterialsNotReplaced().cend())
          {
            // Do not add the current shape material since it should
            // not end up in updateVFs.
            if(name != shape.getMaterial())
              gf_order_by_matnumber.emplace_back(getMaterial(name));
          }
          else
          {
            // The material was in the exclusion list. This means that
            // cannot write to materials that have volume fraction in
            // that zone.
std::cout << "\t" << name << std::endl;
            excludeVFs.emplace_back(getMaterial(name).first);
          }
        }
      }
    }
    // Sort eligible update materials by material number.
    std::sort(gf_order_by_matnumber.begin(), gf_order_by_matnumber.end(),
      [&](const std::pair<mfem::GridFunction *, int> &lhs,
          const std::pair<mfem::GridFunction *, int> &rhs)
    {
      return lhs.second < rhs.second;
    });
std::cout << "*** updateVFs" << std::endl;
    // Append the grid functions in mat number order.
    for(const auto &mat : gf_order_by_matnumber)
    {
std::cout << "\t" << m_vf_material_names[mat.second] << std::endl;
      updateVFs.push_back(mat.first);
    }

    // First time through the shaper, compute VF sums over all materials.
    if(first)
    {
std::cout << "*** Computing initial vf_sums " << std::endl;

      axom::for_all<ExecSpace>(
        dataSize,
        AXOM_LAMBDA(axom::IndexType i) {
          m_vf_sums[i] = 0;
        });
      // Sum each zone's VFs into m_vf_sums.
      int idx = 0;
      for(const auto &gf : m_vf_grid_functions)
      {
std::cout << "*** adding " << m_vf_material_names[idx] << std::endl;

        ArrayView<double,1> matVFView(gf->GetData(), dataSize);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            m_vf_sums[i] += matVFView[i];
          });

gf->Print(std::cout, values_per_line);
idx++;
      }
    }

    // Figure out how much VF in each zone is occupied and immutable.
std::cout << "*** Sum excludeVFs" << std::endl;
    axom::for_all<ExecSpace>(
      dataSize,
      AXOM_LAMBDA(axom::IndexType i) {
        m_vf_occupied[i] = 0.;
      });
    for(auto &gf : excludeVFs)
    {
      ArrayView<double,1> matVFView(gf->GetData(), dataSize);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            m_vf_occupied[i] += matVFView[i];
        });
    }
#if 1
std::cout << "m_vf_occupied={" << std::endl;
for(int i = 0; i < dataSize; i++)
{
    if(i > 0) std::cout << " ";
    std::cout << m_vf_occupied[i];
    if(i > 0 && i % values_per_line == 0)
       std::cout << "\n";
}
std::cout << "}" << std::endl;

std::cout << "*** Writing VFs for shape material" << std::endl;
#endif
    // Compute the volume fractions for the current shape's material.
    ArrayView<double,1> matVFView(matVF.first->GetData(), dataSize);
    axom::for_all<ExecSpace>(
        dataSize,
        AXOM_LAMBDA(axom::IndexType i) {
          // Update this material's VF and m_vf_subtract, which is the
          // amount to subtract from the VF arrays that we need to update.
          double vf = m_overlap_volumes[i] / m_hex_volumes[i];
          if(vf < m_vf_sums[i])
          {
            matVFView[i] = vf;
            m_vf_subtract[i] = 0.;
          }
          else
          {
            double avail = 1. - m_vf_occupied[i];
            double vf_actual = (vf < avail) ? vf : avail;
            matVFView[i] = vf_actual;
            m_vf_subtract[i] = vf_actual;
          }

          // Update the total sum.
          m_vf_sums[i] += matVFView[i] - m_vf_subtract[i];
        });
#if 1
std::cout << "m_vf_sums={" << std::endl;
for(int i = 0; i < dataSize; i++)
{
    if(i > 0) std::cout << " ";
    std::cout << m_vf_sums[i];
    if(i > 0 && i % values_per_line == 0)
       std::cout << "\n";
}
std::cout << "}" << std::endl;
std::cout << "m_vf_subtract={" << std::endl;
for(int i = 0; i < dataSize; i++)
{
    if(i > 0) std::cout << " ";
    std::cout << m_vf_subtract[i];
    if(i > 0 && i % values_per_line == 0)
       std::cout << "\n";
}
std::cout << "}" << std::endl;
std::cout << "matVF={" << std::endl;
matVF.first->Print(std::cout, values_per_line);
std::cout << "}" << std::endl;
#endif
    // Now iterate over updateVFs to subtract off VFs we allocated to the
    // current shape's material.
    for(auto &gf : updateVFs)
    {
      ArrayView<double,1> matVFView(gf->GetData(), dataSize);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            double s = (matVFView[i] < m_vf_subtract[i]) ? matVFView[i] : m_vf_subtract[i];
            matVFView[i] -= s;
            m_vf_subtract[i] -= s;
        });
    }
  }
#endif

  void applyReplacementRules(const klee::Shape& shape) override
  {
    switch(m_execPolicy)
    {
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    case seq:
      applyReplacementRulesImpl<seq_exec>(shape);
      break;
  #if defined(AXOM_USE_OPENMP)
    case omp:
      applyReplacementRulesImpl<omp_exec>(shape);
      break;
  #endif  // AXOM_USE_OPENMP
  #if defined(AXOM_USE_CUDA)
    case cuda:
      applyReplacementRulesImpl<cuda_exec>(shape);
      break;
  #endif  // AXOM_USE_CUDA
  #if defined(AXOM_USE_HIP)
    case hip:
      applyReplacementRulesImpl<hip_exec>(shape);
      break;
  #endif  // AXOM_USE_HIP
#endif    // AXOM_USE_RAJA && AXOM_USE_UMPIRE
    default:
      AXOM_UNUSED_VAR(shape);
      SLIC_ERROR("Unhandled runtime policy case " << m_execPolicy);
      break;
    }
  }
//---------------------------------------------------------------------------
// OLD STUFF
#if 0
  void applyReplacementRules(const klee::Shape& shape) override
  {
    const auto& shapeName = shape.getName();
    const auto& materialName = shape.getMaterial();
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(
        "Applying replacement rules for shape '{}' of material {}",
        shapeName,
        materialName)));

    const std::string vol_frac_("vol_frac_");
    const std::string vol_frac_fmt("vol_frac_{}");

    auto shapeVolFracName = axom::fmt::format("shape_vol_frac_{}", shapeName);
    auto materialVolFracName = axom::fmt::format(vol_frac_fmt, materialName);

    auto* shapeVolFrac = this->getDC()->GetField(shapeVolFracName);
    SLIC_ASSERT(shapeVolFrac != nullptr);

    // Get or create the volume fraction field for this shape's material
    mfem::GridFunction* matVolFrac = nullptr;
    if(this->getDC()->HasField(materialVolFracName))
    {
      matVolFrac = this->getDC()->GetField(materialVolFracName);
    }
    else
    {
      matVolFrac = newVolFracGridFunction();
      this->getDC()->RegisterField(materialVolFracName, matVolFrac);
#if 1
      // Zero out the volume fractions.
      for(int elem = 0; elem < m_num_elements; elem++)
         (*matVolFrac)[elem] = 0.;
#endif
    }

#if 1
#if 1
    // Debugging. Expose our arrays as grid functions.
    static int shaperpass = 0;
#endif
    //-----------------------------------------------------------------------
    // This vector holds the GridFunctions for each material VF object
    // that we need to update (other than that for the current shape,
    // which is referenced using matVolFrac).
    std::vector<mfem::GridFunction *> vfs;

    // Determine the amount that can be written in each element.   
    if(!shape.getMaterialsReplaced().empty())
    {
      // Include only non-shape materials in the replaced list.
      for(const auto &name : shape.getMaterialsReplaced())
      {
        auto vfname = axom::fmt::format(vol_frac_fmt, name);
        if(this->getDC()->HasField(vfname))
        {
          mfem::GridFunction *gf = this->getDC()->GetField(vfname);
          vfs.push_back(gf);
        }
      }

      // These are the VF arrays that have to have something in them for the
      // element to be writeable.
      std::vector<double> writable_vf(m_num_elements, 0.);
      for(auto gf : vfs)
      {
        for(int elem = 0; elem < m_num_elements; elem++)
          writable_vf[elem] += (*gf)[elem];
      }

      // TODO: we could add completely_available here to writable_vf since it is
      //       an additional capacity that could be used for storing the VF.

#if 1
    // Debugging. Expose our arrays as grid functions.
    mfem::GridFunction* ov = newVolFracGridFunction();
    std::string ovname(axom::fmt::format("writable_vf{}", shaperpass));
    shaperpass++;
    this->getDC()->RegisterField(ovname, ov);
    for(int elem = 0; elem < m_num_elements; elem++)
    {
      (*ov)(elem) = writable_vf[elem];
    }
#endif

      for(int elem = 0; elem < m_num_elements; elem++)
      {
        if(writable_vf[elem] > 0.)
        {
          auto VF = m_overlap_volumes[elem] / m_hex_volumes[elem];

          if(VF >= writable_vf[elem] - 1.e-5)
          {
            // There is more VF than writable_vf. We can claim what's there
            // and only that since we can't change other materials in this mode.
            (*matVolFrac)(elem) += writable_vf[elem];
            for(auto gf : vfs)
              (*gf)(elem) = 0.;
          }
          else
          {
            // VF is less than writable_vf, which was the amount that other
            // materials that we're replacing occupy. We need to make
            // all of the materials fit in writable_vf.
            double scale = (writable_vf[elem] - VF) / writable_vf[elem];
            for(auto gf : vfs)
              (*gf)(elem) *= scale;
            (*matVolFrac)(elem) = VF;
          }
        }
      }
    }
    else if(!shape.getMaterialsNotReplaced().empty())
    {
      // Include all materials that are NOT in replaced list.
      std::vector<mfem::GridFunction *> forbidden;
      std::vector<std::string> forbidden_names;
      for(auto it : this->getDC()->GetFieldMap())
      {
        if(it.first.find(vol_frac_) == 0)
        {
          std::string name(it.first.substr(vol_frac_.size()));
          // See if the material name is in the list of materials NOT replaced.
          auto it2 = std::find(shape.getMaterialsNotReplaced().cbegin(),
                               shape.getMaterialsNotReplaced().cend(),
                               name);
          if(it2 != shape.getMaterialsNotReplaced().cend())
          {
            forbidden.push_back(it.second);
            forbidden_names.push_back(it.first);
          }
          else
            vfs.push_back(it.second);
        }
      }
      // NOTE: the shape material is included in vfs so we can include it
      //       in writable_vf.

      // Assume we can write everywhere. Then we subtract off the forbidden areas.
      std::vector<double> writable_vf(m_num_elements, 1.);
      for(size_t i = 0; i < forbidden.size(); i++)
      {
        auto gf = forbidden[i];
        std::cout << "!!!!!!!!!!!!!! " << forbidden_names[i] << " is FORBIDDEN!" << std::endl;
        for(int elem = 0; elem < m_num_elements; elem++)
          writable_vf[elem] -= (*gf)[elem];
      }
      // Make sure we get zeroes if things are really close to zero.
      for(int elem = 0; elem < m_num_elements; elem++)
      {
          writable_vf[elem] = abs(writable_vf[elem]);
          if(writable_vf[elem] < 1.e-10)
             writable_vf[elem] = 0.;
      }

#if 1
    // Debugging. Expose our arrays as grid functions.
    mfem::GridFunction* ov = newVolFracGridFunction();
    std::string ovname(axom::fmt::format("writable_vf{}", shaperpass));
    shaperpass++;
    this->getDC()->RegisterField(ovname, ov);
    for(int elem = 0; elem < m_num_elements; elem++)
    {
      (*ov)(elem) = writable_vf[elem];
    }
#endif
    int skips = 0;
      for(int elem = 0; elem < m_num_elements; elem++)
      {
        if(writable_vf[elem] > 0.0)
        {
          auto VF = m_overlap_volumes[elem] / m_hex_volumes[elem];
          if(VF >= writable_vf[elem])
          {
            // There is more VF than writable_vf. We can claim what's there
            // and only that since we can't change other materials in this mode.
            for(auto gf : vfs)
              (*gf)(elem) = 0.;
            (*matVolFrac)(elem) = writable_vf[elem];
          }
          else
          {
            // VF is less than writable_vf, which was the amount that other
            // materials that we're replacing occupy. We need to make
            // all of the materials fit in writable_vf.
            double scale = (writable_vf[elem] - VF) / writable_vf[elem];
            for(auto gf : vfs)
              (*gf)(elem) *= scale;
            (*matVolFrac)(elem) = VF;
          }
        }
        else
        {
           skips++;
           (*matVolFrac)(elem) = 0.;
        }
      }

      std::cout << "!!!!!!!!!!!!!! Skipped " << skips << " elements when writing material " << shape.getMaterial() << std::endl;
    }
    else
    {
      // Normal case. Every non-shape material is allowed to be replaced.
      for(auto it : this->getDC()->GetFieldMap())
      {
        if(it.first.find(vol_frac_) == 0)
        {
          std::string name(it.first.substr(vol_frac_.size()));
          if(name == shape.getMaterial())
            continue;
          vfs.push_back(it.second);
        }
      }

      // For each cell, figure out how much VF is still completely free
      // (unallocated to any material).
      std::vector<double> writable_vf(m_num_elements, 1.);
      for(auto it : this->getDC()->GetFieldMap())
      {
        if(it.first.find(vol_frac_) == 0)
        {
          const mfem::GridFunction *gf = it.second;          
          for(int elem = 0; elem < m_num_elements; elem++)
          {
            writable_vf[elem] -= (*gf)[elem];
          }
        }
      }
      // Clamp any negative epsilon values.
      for(int elem = 0; elem < m_num_elements; elem++)
      {
        writable_vf[elem] = std::max(writable_vf[elem], 0.);
      }

      for(int elem = 0; elem < m_num_elements; elem++)
      {
          auto VF = m_overlap_volumes[elem] / m_hex_volumes[elem];
          constexpr double ONE_TOL = 1. - 1.e-10;
          if(VF >= ONE_TOL)
          {
            // This material gets it all.
            (*matVolFrac)(elem) = 1.;
            for(auto gf : vfs)
              (*gf)(elem) = 0.;
          }
          else if(VF <= writable_vf[elem])
          {
            // This material can claim all of its VF because the element
            // still has enough void to fit it.
            (*matVolFrac)(elem) += VF;
          }
          else
          {
            // VF is less than writable_vf, which was the amount that other
            // materials that we're replacing occupy. We need to make
            // all of the materials fit in writable_vf.
            double other_vf = 1. - writable_vf[elem];
            double scale = (other_vf - VF) / other_vf;

            // If writable_vf[elem] is 0. then other_vf will be 1.
            // scale = (1 - VF)

            for(auto gf : vfs)
              (*gf)(elem) *= scale;
            (*matVolFrac)(elem) = VF;
          }
      }
    }

#endif

#if 0
#if 1
//--------------------------------------------------------------------------------------
    // For each cell, figure out how much VF is still completely free
    // (unallocated to any material).
    std::vector<double> completely_free(m_num_elements, 1.);
    for(auto it : this->getDC()->GetFieldMap())
    {
      if(it.first.find(vol_frac_) == 0)
      {
        const mfem::GridFunction *gf = it.second;          
        for(int elem = 0; elem < m_num_elements; elem++)
        {
          completely_free[elem] -= (*gf)[elem];
        }
      }
    }
    // Clamp any negative epsilon values.
    for(int elem = 0; elem < m_num_elements; elem++)
    {
      completely_free[elem] = std::max(completely_free[elem], 0.);
    }

    // This vector holds the GridFunctions for each material VF object
    // that we need to update (other than that for the current shape,
    // which is referenced using matVolFrac).
    std::vector<mfem::GridFunction *> vfs;

    // Populate vfs with the pointers to the grid functions that can be
    // updated when we modify VFs. This does not include the grid function
    // for the current shape, matVolFrac.
    if(!shape.getMaterialsReplaced().empty())
    {
      // Include only non-shape materials in the replaced list.
      for(const auto &name : shape.getMaterialsReplaced())
      {
        if(name == shape.getMaterial())
          continue;
        auto vfname = axom::fmt::format(vol_frac_fmt, name);
        if(this->getDC()->HasField(vfname))
        {
          mfem::GridFunction *gf = this->getDC()->GetField(vfname);
          vfs.push_back(gf);
        }
      }
    }
    else if(!shape.getMaterialsNotReplaced().empty())
    {
      // We include all non-shape materials unless they are in the NOT
      // replaced list.
      for(auto it : this->getDC()->GetFieldMap())
      {
        if(it.first.find(vol_frac_) == 0)
        {
          std::string name(it.first.substr(vol_frac_.size()));
          if(name == shape.getMaterial())
            continue;
          // See if the material name is in the list of materials NOT replaced.
          auto it2 = std::find(shape.getMaterialsNotReplaced().cbegin(),
                               shape.getMaterialsNotReplaced().cend(),
                               name);
          if(it2 == shape.getMaterialsNotReplaced().cend())
            vfs.push_back(it.second);
        }
      }
    }
    else
    {
      // Normal case. Every non-shape material is allowed to be replaced.
      for(auto it : this->getDC()->GetFieldMap())
      {
        if(it.first.find(vol_frac_) == 0)
        {
          std::string name(it.first.substr(vol_frac_.size()));
          if(name == shape.getMaterial())
            continue;
          vfs.push_back(it.second);
        }
      }
    }

    // Print some debugging info about which replacements will happen.
    for(auto it : this->getDC()->GetFieldMap())
    {
      if(it.first.find(vol_frac_) == 0)
      {
        std::string name(it.first.substr(vol_frac_.size()));
        bool rep = std::find(vfs.begin(), vfs.end(), it.second) != vfs.cend();
        SLIC_DEBUG(axom::fmt::format(
          "Should we replace material '{}' with shape '{}' of material '{}'? {}",
          name,
          shape.getName(),
          shape.getMaterial(),
          rep ? "yes" : "no"));
      }
    }

    // For each cell, figure out how much VF can be overwritten if
    // we need the space. This is the sum of all VFs that we are
    // allowed to write into.
    std::vector<double> overwriteable(m_num_elements, 0.);
    for(int elem = 0; elem < m_num_elements; elem++)
    {
      overwriteable[elem] += (*matVolFrac)(elem);
    }
    for(const auto gf : vfs)
    {
      for(int elem = 0; elem < m_num_elements; elem++)
      {
        overwriteable[elem] += (*gf)(elem);
      }
    }

#if 1
    // Debugging. Expose our arrays as grid functions.
    static int shaperpass = 0;
    mfem::GridFunction* cf = newVolFracGridFunction();
    mfem::GridFunction* ov = newVolFracGridFunction();
    std::string cfname(axom::fmt::format("cf{}", shaperpass));
    std::string ovname(axom::fmt::format("ov{}", shaperpass));
    shaperpass++;
    this->getDC()->RegisterField(cfname, cf);
    this->getDC()->RegisterField(ovname, ov);
    for(int elem = 0; elem < m_num_elements; elem++)
    {
      (*cf)(elem) = completely_free[elem];
      (*ov)(elem) = overwriteable[elem];
    }
#endif

    // Now update the volume fractions.
    for(int elem = 0; elem < m_num_elements; elem++)
    {
      auto VF = m_overlap_volumes[elem] / m_hex_volumes[elem];
      if(VF <= completely_free[elem])
      {
        // There is enough free room for this VF. No need to update other VFs.
        (*matVolFrac)(elem) += VF;
      }
      else
      {
        double avail = completely_free[elem] + overwriteable[elem];
        constexpr double ONE_TOL = 1. - 1.e-10;
        if(VF >= ONE_TOL)
        {
          // VF takes the whole element.
          for(auto &gf : vfs)
             (*gf)(elem) = 0.;
          (*matVolFrac)(elem) = VF;
        }
        else if(VF >= avail)
        {
          // VF can completely overtake the space in the available VFs.
          // This material gets it all and we zero out the others.
          for(auto &gf : vfs)
             (*gf)(elem) = 0.;
          (*matVolFrac)(elem) = avail;
        }
        else
        {
          // We can fit the VF in avail but it does not cover all of avail.
          // So, we?ll scale the other replacement materials
          double scale = avail / (overwriteable[elem] + VF);
          (*matVolFrac)(elem) = scale * VF;
          // Scale the replaced materials.
          for(auto gf : vfs)
             (*gf)(elem) *= scale;
        }
      }
    }

    // Update the shape volume fractions.
    for(int i = 0; i < m_num_elements; i++)
    {
      (*shapeVolFrac)(i) = m_overlap_volumes[i] / m_hex_volumes[i];
    }
#else
    // update material volume fractions
    for(int i = 0; i < m_num_elements; i++)
    {
      (*matVolFrac)(i) = m_overlap_volumes[i] / m_hex_volumes[i];
    }

    /// Implementation here -- update material volume fractions based on replacement rules
    // Note: we're not yet updating the shape volume fractions
    AXOM_UNUSED_VAR(shapeVolFrac);
#endif
#endif
  }
#endif

  void finalizeShapeQuery() override
  {
    // Implementation here -- destroy BVH tree and other shape-based data structures
    delete m_surfaceMesh;
    axom::deallocate(m_hex_volumes);
    axom::deallocate(m_overlap_volumes);

    m_surfaceMesh = nullptr;
  }

  //@}

public:
  // Prepares the shaping query, based on the policy member set
  // (default is sequential)
  void prepareShapeQuery(klee::Dimensions shapeDimension,
                         const klee::Shape& shape) override
  {
    switch(m_execPolicy)
    {
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    case seq:
      prepareShapeQueryImpl<seq_exec>(shapeDimension, shape);
      break;
  #if defined(AXOM_USE_OPENMP)
    case omp:
      prepareShapeQueryImpl<omp_exec>(shapeDimension, shape);
      break;
  #endif  // AXOM_USE_OPENMP
  #if defined(AXOM_USE_CUDA)
    case cuda:
      prepareShapeQueryImpl<cuda_exec>(shapeDimension, shape);
      break;
  #endif  // AXOM_USE_CUDA
  #if defined(AXOM_USE_HIP)
    case hip:
      prepareShapeQueryImpl<hip_exec>(shapeDimension, shape);
      break;
  #endif  // AXOM_USE_HIP
#endif    // AXOM_USE_RAJA && AXOM_USE_UMPIRE
    default:
      AXOM_UNUSED_VAR(shapeDimension);
      AXOM_UNUSED_VAR(shape);
      SLIC_ERROR("Unhandled runtime policy case " << m_execPolicy);
      break;
    }
  }

  // Runs the shaping query, based on the policy member set
  // (default is sequential)
  void runShapeQuery(const klee::Shape& shape) override
  {
    switch(m_execPolicy)
    {
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    case seq:
      runShapeQueryImpl<seq_exec>(shape);
      break;
  #if defined(AXOM_USE_OPENMP)
    case omp:
      runShapeQueryImpl<omp_exec>(shape);
      break;
  #endif  // AXOM_USE_OPENMP
  #if defined(AXOM_USE_CUDA)
    case cuda:
      runShapeQueryImpl<cuda_exec>(shape);
      break;
  #endif  // AXOM_USE_CUDA
  #if defined(AXOM_USE_HIP)
    case hip:
      runShapeQueryImpl<hip_exec>(shape);
      break;
  #endif  // AXOM_USE_HIP
#endif    // AXOM_USE_RAJA && AXOM_USE_UMPIRE
    default:
      AXOM_UNUSED_VAR(shape);
      SLIC_ERROR("Unhandled runtime policy case " << m_execPolicy);
      break;
    }
  }

  void adjustVolumeFractions() override
  {
    // Implementation here -- not sure if this will require anything for intersection-based shaping
  }

private:
  /// Create and return a new volume fraction grid function for the current mesh
  mfem::GridFunction* newVolFracGridFunction()
  {
    mfem::Mesh* mesh = getDC()->GetMesh();
    SLIC_ASSERT(mesh != nullptr);

    const int vfOrder = 0;
    const int dim = mesh->Dimension();
    mfem::L2_FECollection* coll =
      new mfem::L2_FECollection(vfOrder, dim, mfem::BasisType::Positive);
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, coll);
    mfem::GridFunction* volFrac = new mfem::GridFunction(fes);
    volFrac->MakeOwner(coll);

    return volFrac;
  }

private:
  ExecPolicy m_execPolicy {seq};
  int m_level {7};
  int m_num_elements {0};
  double* m_hex_volumes {nullptr};
  double* m_overlap_volumes {nullptr};
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  double m_vertexWeldThreshold {1.e-10};
  int m_octcount {0};
  OctahedronType* m_octs {nullptr};
  BoundingBoxType* m_aabbs {nullptr};
  PolyhedronType* m_hexes {nullptr};
  BoundingBoxType* m_hex_bbs {nullptr};

  int m_replacement_allocator{-1};
  double* m_vf_sums {nullptr};
  double* m_vf_subtract {nullptr};
  double* m_vf_occupied {nullptr};
  std::vector<mfem::GridFunction *> m_vf_grid_functions;
  std::vector<std::string> m_vf_material_names;

#endif
  // What do I need here?
  // Probably size of stuff
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_INTERSECTION_SHAPER__HPP_
