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
    slic::flushStreams();
    umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

    // Save current/default allocator
    const int current_allocator = axom::getDefaultAllocatorID();

    // Determine new allocator (for CUDA/HIP policy, set to Unified)
    // Set new default to device
    axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

    const auto& shapeName = shape.getName();
    AXOM_UNUSED_VAR(shapeDimension);
    AXOM_UNUSED_VAR(shapeName);

    SLIC_INFO(axom::fmt::format("Current shape is {}", shapeName));
slic::flushStreams();
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
slic::flushStreams();
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
slic::flushStreams();
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

    SLIC_INFO("!!!!!!!!!!!!!!!!!!!!!!!!!");
    slic::flushStreams();
    // Generate the Octahedra
    const bool disc_status = axom::quest::discretize<ExecSpace>(polyline,
                                                                polyline_size,
                                                                m_level,
                                                                m_octs,
                                                                m_octcount);

    // Oddities required by hip
    OctahedronType* local_octs = m_octs;

    AXOM_UNUSED_VAR(disc_status);  // silence warnings in release configs
    SLIC_ASSERT_MSG(
      disc_status,
      "Discretization of contour has failed. Check that contour is valid");

    SLIC_INFO(
      axom::fmt::format("Contour has been discretized into {} octahedra ",
                        m_octcount));
slic::flushStreams();
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
      slic::flushStreams();

      // Print out the total volume of all the octahedra
      using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
      RAJA::ReduceSum<REDUCE_POL, double> total_oct_vol(0.0);
      axom::for_all<ExecSpace>(
        m_octcount,
        AXOM_LAMBDA(axom::IndexType i) {
          // Convert Octahedron into Polyhedrom
          PolyhedronType octPoly;
          double octVolume;

          // OK
          // octPoly.addVertex(Point3D{});

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
slic::flushStreams();

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
slic::flushStreams();

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

    umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();

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
slic::flushStreams();
    // Generate the BVH tree over the octahedra
    // Access-aligned bounding boxes
    m_aabbs = axom::allocate<BoundingBoxType>(m_octcount);
    // BoundingBoxType * local_aabbs = axom::allocate<BoundingBoxType>(m_octcount);

    // m_aabbs[0] = BoundingBoxType{};

    // OK
    // BoundingBoxType* lebox = axom::allocate<BoundingBoxType>(m_octcount);

    // Oddities required by hip
    OctahedronType* local_octs = m_octs;
    BoundingBoxType * local_aabbs = m_aabbs;

    SLIC_INFO("AAAAAAAA");
    // Get the bounding boxes for the Octahedrons
    axom::for_all<ExecSpace>(
      m_octcount,
      AXOM_LAMBDA(axom::IndexType i) {
        // m_aabbs[i] = primal::compute_bounding_box<double, 3>(m_octs[i]);
        local_aabbs[i] = primal::compute_bounding_box<double, 3>(local_octs[i]);

        // OK.
        // local_octs[i] = local_octs[i];

        // OK
        // lebox[i].isValid();

        // This works
        // BoundingBoxType doggo;
        // OctahedronType oxa;
        // doggo = primal::compute_bounding_box<double, 3>(oxa);
      });

    SLIC_INFO("BBBBBBBBBBBB");

    // Insert Octahedra Bounding Boxes into BVH.
    //bvh.setAllocatorID(poolID);
    spin::BVH<3, ExecSpace, double> bvh;
    bvh.initialize(m_aabbs, m_octcount);

        SLIC_INFO("CCCCCCCCCC");

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

SLIC_INFO("DDDDDDDDDDDDDDDD");
    // Create and register a scalar field for this shape's volume fractions
    // The Degrees of Freedom will be in correspondence with the elements
    auto* volFrac = this->newVolFracGridFunction();
    auto volFracName = axom::fmt::format("shape_vol_frac_{}", shape.getName());
    this->getDC()->RegisterField(volFracName, volFrac);

SLIC_INFO("EEEEEEEEEEEEEEE");
    // Initialize hexahedral elements
    m_hexes = axom::allocate<PolyhedronType>(NE);
    m_hex_bbs = axom::allocate<BoundingBoxType>(NE);

    // Oddities required by hip
    PolyhedronType * local_hexes = m_hexes;
    BoundingBoxType * local_hex_bbs = m_hex_bbs; 

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

SLIC_INFO("FFFFFFFFFFFFF");
slic::flushStreams();

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
          if(axom::utilities::isNearlyEqual(local_hexes[i][j][0], 0.0, ZERO_THRESHOLD))
          {
            local_hexes[i][j][0] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_hexes[i][j][1], 0.0, ZERO_THRESHOLD))
          {
            local_hexes[i][j][1] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_hexes[i][j][2], 0.0, ZERO_THRESHOLD))
          {
            local_hexes[i][j][2] = 0.0;
          }
        }

        // Get bounding box for hexahedral element
        local_hex_bbs[i] = primal::compute_bounding_box<double, 3>(local_hexes[i]);
      });  // end of loop to initialize hexahedral elements and bounding boxes

SLIC_INFO("GGGGGGGGGGGGG");
    // Deallocatem no longer needed
    axom::deallocate(vertCoords);

    // Set octahedra components to zero if within threshold
    axom::for_all<ExecSpace>(
      m_octcount,
      AXOM_LAMBDA(axom::IndexType i) {
        for(int j = 0; j < OctahedronType::NUM_OCT_VERTS; j++)
        {
          if(axom::utilities::isNearlyEqual(local_octs[i][j][0], 0.0, ZERO_THRESHOLD))
          {
            local_octs[i][j][0] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_octs[i][j][1], 0.0, ZERO_THRESHOLD))
          {
            local_octs[i][j][1] = 0.0;
          }

          if(axom::utilities::isNearlyEqual(local_octs[i][j][2], 0.0, ZERO_THRESHOLD))
          {
            local_octs[i][j][2] = 0.0;
          }
        }
      });

SLIC_INFO("HHHHHHHHHHHH");
slic::flushStreams();

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

    // Oddities required by hip
    double * local_overlap_volumes = m_overlap_volumes;
    double * local_hex_volumes = m_hex_volumes;

SLIC_INFO("IIIIIIIIIIIIIII");
slic::flushStreams();

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
          PolyhedronType poly = primal::clip(local_octs[octIndex], tets[tetIndex]);

          // Poly is valid
          if(poly.numVertices() >= 4)
          {
            double clip_volume = poly.volume();
            // Flip sign if negative
            if(clip_volume < 0)
            {
              clip_volume = -clip_volume;
            }
            RAJA::atomicAdd<ATOMIC_POL>(local_overlap_volumes + index, clip_volume);
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
slic::flushStreams();


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

    auto shapeVolFracName = axom::fmt::format("shape_vol_frac_{}", shapeName);
    auto materialVolFracName = axom::fmt::format("vol_frac_{}", materialName);

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
    }

    // update material volume fractions
    for(int i = 0; i < m_num_elements; i++)
    {
      (*matVolFrac)(i) = m_overlap_volumes[i] / m_hex_volumes[i];
    }

    /// Implementation here -- update material volume fractions based on replacement rules
    // Note: we're not yet updating the shape volume fractions
    AXOM_UNUSED_VAR(shapeVolFrac);
  }

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
    SLIC_INFO("INSIDE prepareShapeQuery");
    slic::flushStreams();

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
        SLIC_INFO("INSIDE prepareShapeQuery HIP HIP HIP");
    slic::flushStreams();
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
#endif
  // What do I need here?
  // Probably size of stuff
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_INTERSECTION_SHAPER__HPP_
