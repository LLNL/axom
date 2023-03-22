// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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
/*!
 * \class GridFunctionView
 *
 * \brief Provides a view over an MFEM grid function. MFEM grid functions are
 *        assumed to live in host memory. This class performs data movement
 *        needed to access the grid function data within a GPU device lambda. This
 *        view is limited in scope, though could be expanded in the future.
 *
 * \tparam ExecSpace The execution space where the grid function data will
 *                   be accessed.
 */
template <typename ExecSpace>
class GridFunctionView
{
public:
  /*!
   * \brief Host constructor that accepts the grid function.
   *
   * \param gf The grid function that will be accessed/modified by the view.
   * \param _needResult Whether the data needs to be brought back to the host
   *                    from the device.
   */
  AXOM_HOST GridFunctionView(mfem::GridFunction* gf, bool _needResult = true)
  {
    initialize(gf->GetData(), gf->Size(), _needResult);
  }

  /*!
   * \brief Copy constructor, which is called to make a copy of the host
   *        object so it is accessible inside a RAJA kernel. Any data movement
   *        happened in the host constructor. This version sets hostData to
   *        nullptr so we know not to clean up in the destructor.
   */
  AXOM_HOST_DEVICE GridFunctionView(const GridFunctionView& obj)
    : m_hostData(nullptr)
    , m_deviceData(obj.m_deviceData)
    , m_numElements(obj.m_numElements)
    , m_needResult(obj.m_needResult)
  { }

  /*!
   * \brief Destructor. On the host, this method may move data from the 
            device and deallocate device storage.
   */
  AXOM_HOST_DEVICE ~GridFunctionView() { finalize(); }

  /*!
   * \brief Indexing operator for accessing the data.
   *
   * \param i The index at which to access the data.
   *
   * \return A reference to the data at index i.
   */
  AXOM_HOST_DEVICE double& operator[](int i) { return m_deviceData[i]; }
  // non-const return on purpose.
  AXOM_HOST_DEVICE double& operator[](int i) const { return m_deviceData[i]; }

private:
  /*!
   * \brief Initializes members using data from the grid function. This method
   *        is called on the host.
   *
   * \param hostPtr The grid function data pointer on the host.
   * \param nElem   The grid function size.
   * \param _needResult Whether any data are copied from device.
   */
  AXOM_HOST void initialize(double* hostPtr, int nElem, bool _needResult)
  {
    m_hostData = m_deviceData = hostPtr;
    m_numElements = nElem;
    m_needResult = _needResult;
  }

  /*!
   * \brief Helps during destruction.
   */
  AXOM_HOST_DEVICE void finalize() { m_deviceData = nullptr; }

#if defined(AXOM_USE_CUDA) || defined(AXOM_USE_HIP)
  /*!
   * \brief Initializes members using data from the grid function. This method
   *        is called on the host and it copies data to the device.
   *
   * \param hostPtr The grid function data pointer on the host.
   * \param nElem   The grid function size.
   * \param _needResult Whether any data are copied from device.
   */
  AXOM_HOST void initializeDevice(double* hostPtr, int nElem, bool _needResult)
  {
    m_hostData = hostPtr;
    m_numElements = nElem;
    m_needResult = _needResult;
    int execSpaceAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto dataSize = sizeof(double) * m_numElements;
    m_deviceData = axom::allocate<double>(dataSize, execSpaceAllocatorID);
    axom::copy(m_deviceData, m_hostData, dataSize);
  }

  /*!
   * \brief Helps during destruction. On the host, it copies device data back
   *        into the grid function on the host.
   */
  AXOM_HOST_DEVICE void finalizeDevice()
  {
  #ifndef AXOM_DEVICE_CODE
    // Only the host will do this work.
    if(m_hostData != nullptr)
    {
      if(m_needResult)
      {
        auto dataSize = sizeof(double) * m_numElements;
        axom::copy(m_hostData, m_deviceData, dataSize);
      }
      axom::deallocate(m_deviceData);
      m_deviceData = nullptr;
    }
  #endif
  }
#endif

private:
  double* m_hostData {nullptr};
  double* m_deviceData {nullptr};
  int m_numElements {0};
  bool m_needResult {false};
};

#if defined(AXOM_USE_CUDA)
/*!
 * \brief CUDA specialization that calls initializeDevice to copy data
 *        from the host to the device.
 *
 * \param hostPtr The grid function data pointer on the host.
 * \param nElem   The grid function size.
 * \param _needResult Whether any data are copied from device.
 */
template <>
AXOM_HOST inline void GridFunctionView<cuda_exec>::initialize(double* hostPtr,
                                                              int nElem,
                                                              bool _needResult)
{
  initializeDevice(hostPtr, nElem, _needResult);
}

/*!
 * \brief CUDA specialization that may copy data back from the device
 *        and deallocate any associated device data.
 */
template <>
AXOM_HOST_DEVICE inline void GridFunctionView<cuda_exec>::finalize()
{
  finalizeDevice();
}
#endif
#if defined(AXOM_USE_HIP)
/*!
 * \brief HIP specialization that calls initializeDevice to copy data
 *        from the host to the device.
 *
 * \param hostPtr The grid function data pointer on the host.
 * \param nElem   The grid function size.
 * \param _needResult Whether any data are copied from device.
 */
template <>
AXOM_HOST inline void GridFunctionView<hip_exec>::initialize(double* hostPtr,
                                                             int nElem,
                                                             bool _needResult)
{
  initializeDevice(hostPtr, nElem, _needResult);
}

/*!
 * \brief HIP specialization that may copy data back from the device
 *        and deallocate any associated device data.
 */
template <>
AXOM_HOST_DEVICE inline void GridFunctionView<hip_exec>::finalize()
{
  finalizeDevice();
}
#endif

//---------------------------------------------------------------------------
class IntersectionShaper : public Shaper
{
public:
  using BoundingBoxType = primal::BoundingBox<double, 3>;
  using HexahedronType = primal::Hexahedron<double, 3>;
  using OctahedronType = primal::Octahedron<double, 3>;
  using PolyhedronType = primal::Polyhedron<double, 3>;
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
    m_hexes = axom::allocate<HexahedronType>(NE);
    m_hex_bbs = axom::allocate<BoundingBoxType>(NE);

    // Oddities required by hip to avoid capturing `this`
    HexahedronType* local_hexes = m_hexes;
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
        local_hexes[i] = HexahedronType();
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIndex = (i * NUM_VERTS_PER_HEX * NUM_COMPS_PER_VERT) +
            j * NUM_COMPS_PER_VERT;
          local_hexes[i][j] = Point3D({vertCoords[vertIndex],
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

    using TetHexArray = axom::StackArray<TetrahedronType, NUM_TETS_PER_HEX>;

    AXOM_PERF_MARK_SECTION("init_tets",
                           axom::for_all<ExecSpace>(
                             NE,
                             AXOM_LAMBDA(axom::IndexType i) {
                               TetHexArray cur_tets;
                               local_hexes[i].triangulate(cur_tets);

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
                             NE,
                             AXOM_LAMBDA(axom::IndexType i) {
                               local_hex_volumes[i] = local_hexes[i].volume();
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
  /*!
   * \brief Turn a material name into a grid function name.
   *
   * \param materialName The name of the material.
   *
   * \return The name of the material's grid function.
   */
  std::string materialNameToFieldName(const std::string& materialName) const
  {
    const std::string vol_frac_fmt("vol_frac_{}");
    auto name = axom::fmt::format(vol_frac_fmt, materialName);
    return name;
  }

  /*!
   * \brief Turn a grid function name into a material name.
   *
   * \param fieldName The name of the grid function.
   *
   * \return The name of the material material.
   */
  std::string fieldNameToMaterialName(const std::string& fieldName) const
  {
    const std::string vol_frac_("vol_frac_");
    std::string name;
    if(fieldName.find(vol_frac_) == 0)
      name = fieldName.substr(vol_frac_.size());
    return name;
  }

  /*!
   * \brief Gets the grid function and material number for a material name.
   *
   * \param materialName The name of the material.
   *
   * \return A pair containing the associated grid function and material
   *         number (its order in the list).
   */
  std::pair<mfem::GridFunction*, int> getMaterial(const std::string& materialName)
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
    }
    else
    {
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

  /*!
   * \brief Scans the grid functions in the data collection and creates
   *        a material entry for any that do not already exist. We maintain
   *        our own vectors because we assume that the order of materials
   *        does not change. The mfem::DataCollection uses a map internally
   *        so if we add materials, it could change the traversal order.
   *
   */
  void populateMaterials()
  {
    std::vector<std::string> materialNames;
    for(auto it : this->getDC()->GetFieldMap())
    {
      std::string materialName = fieldNameToMaterialName(it.first);
      if(!materialName.empty()) materialNames.emplace_back(materialName);
    }
    // Add any of these existing fields to this class' bookkeeping.
    for(const auto& materialName : materialNames)
    {
      (void)getMaterial(materialName);
    }
  }

  // Switch back to public. This is done here because the CUDA compiler
  // does not like the following template functions to be private.
public:
  /*!
   * \brief Make a new grid function that contains all of the free space not
   *        occupied by existing materials.
   *
   * \note We currently leave completely_free in the data collection,
   *       even after the shaper has executed.
   *
   * \tparam ExecSpace The execution space where the data are computed.
   *
   * \return The grid function that represents the amount of completely
   *         free space in each zone.
   */
  template <typename ExecSpace>
  mfem::GridFunction* getCompletelyFree()
  {
    const std::string fieldName("completely_free");
    mfem::GridFunction* cfgf = nullptr;
    if(this->getDC()->HasField(fieldName))
      cfgf = this->getDC()->GetField(fieldName);
    else
    {
      // Make the new grid function.
      cfgf = newVolFracGridFunction();
      this->getDC()->RegisterField(fieldName, cfgf);

      AXOM_PERF_MARK_SECTION("compute_completely_free", {
        int dataSize = cfgf->Size();
        GridFunctionView<ExecSpace> cfView(cfgf);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) { cfView[i] = 1.; });

        // Iterate over all materials and subtract off their VFs from cfgf.
        for(auto& gf : m_vf_grid_functions)
        {
          GridFunctionView<ExecSpace> matVFView(gf, false);
          axom::for_all<ExecSpace>(
            dataSize,
            AXOM_LAMBDA(axom::IndexType i) {
              cfView[i] -= matVFView[i];
              cfView[i] = (cfView[i] < 0.) ? 0. : cfView[i];
            });
        }
      });
    }
    return cfgf;
  }

  /*!
   * \brief Set the volume fractions for the current shape into the grid
   *        function for the material and adjust any other material volume
   *        fraction grid functions.
   *
   * \tparam ExecSpace The execution space where the data are computed.
   * \param shape The shape whose volume fractions are being stored.
   *
   * The replacement rules operate on the current shape's material as well as the
   * rest of the materials since they need to be updated to sum to 1. When adding
   * a volume fraction to a zone, the code adds the allowed amount to the zone.
   * In most cases, this is just the computed material volume, though it can be
   * restricted by amounts of other materials if replacing materials. Once a
   * volume fraction has been added, a list of update materials is iterated and
   * a corresponding amount is subtracted from them. This update list consists of
   * "completely free", all non-shape  materials, and finally the shape material
   * itself. The shape material is added at the end in case the initial volume
   * fraction addition exceeded 1.
   *
   * The "completely free" material (CF) represents any volume fraction in the
   * zones that has not been assigned to any material. We subtract from CF first
   * as part of the update process so we do not have to subtract from real materials
   * in the event that a zone is not full.
   *
   * The list of update materials depends on the shape's material replacement rule.
   *
   * Example:
   *                        update mats
   *                        |---->
   *
   *                  mat1  |  CF     mat0
   *                  ---------------------
   * add 0.4 mat1  -> | 0.4 | 0.2  | 0.4  |
   *                  ---------------------
   *
   *                  mat1     CF     mat0     subtract
   *                  ---------------------    -------
   *                  | 0.8 | 0.2  | 0.4  |    | 0.4 |  (added 0.4 to mat1)
   *                  ---------------------    -------
   *
   *                  mat1     CF     mat0     subtract
   *                  ---------------------    -------
   *                  | 0.8 | 0.0  | 0.4  |    | 0.2 |  (subtracted 0.2 from CF)
   *                  ---------------------    -------
   *
   *                  mat1     CF     mat0     subtract
   *                  ---------------------    -------
   *                  | 0.8 | 0.0  | 0.2  |    | 0.0 |  (subtracted 0.2 from mat0)
   *                  ---------------------    -------
   */
  template <typename ExecSpace>
  void applyReplacementRulesImpl(const klee::Shape& shape)
  {
    // Make sure the material lists are up to date.
    populateMaterials();

    // Get this shape's material, creating the GridFunction if needed.
    auto matVF = getMaterial(shape.getMaterial());
    int dataSize = matVF.first->Size();

    // Get this shape's array.
    auto shapeVolFracName =
      axom::fmt::format("shape_vol_frac_{}", shape.getName());
    auto* shapeVolFrac = this->getDC()->GetField(shapeVolFracName);
    SLIC_ASSERT(shapeVolFrac != nullptr);

    // Allocate some memory for the replacement rule data arrays.
    int execSpaceAllocatorID = axom::execution_space<ExecSpace>::allocatorID();
    Array<double> vf_subtract_array(dataSize, dataSize, execSpaceAllocatorID);
    Array<double> vf_writable_array(dataSize, dataSize, execSpaceAllocatorID);
    ArrayView<double> vf_subtract(vf_subtract_array);
    ArrayView<double> vf_writable(vf_writable_array);

    // Determine which grid functions need to be considered for VF updates.
    std::vector<std::pair<mfem::GridFunction*, int>> gf_order_by_matnumber;
    std::vector<mfem::GridFunction*> updateVFs, excludeVFs;
    if(!shape.getMaterialsReplaced().empty())
    {
      // Include materials replaced in updateVFs.
      for(const auto& name : shape.getMaterialsReplaced())
      {
        gf_order_by_matnumber.emplace_back(getMaterial(name));
      }
    }
    else
    {
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
            excludeVFs.emplace_back(getMaterial(name).first);
          }
        }
      }
    }
    // Sort eligible update materials by material number.
    std::sort(gf_order_by_matnumber.begin(),
              gf_order_by_matnumber.end(),
              [&](const std::pair<mfem::GridFunction*, int>& lhs,
                  const std::pair<mfem::GridFunction*, int>& rhs) {
                return lhs.second < rhs.second;
              });

    // Append the completely free grid function to the materials we update
    // Add it first so it is the highest priority material. This helps us
    // account for how much we need to deduct from the real materials during
    // subtraction.
    updateVFs.push_back(getCompletelyFree<ExecSpace>());

    // Append the grid functions in mat number order.
    for(const auto& mat : gf_order_by_matnumber)
    {
      updateVFs.push_back(mat.first);
    }

    // We put the current shape's material at the end of the update list,
    // in the case that adding to its VF initially exceeds 1.
    updateVFs.push_back(getMaterial(shape.getMaterial()).first);

    // Figure out how much VF in each zone can be written.
    if(!shape.getMaterialsReplaced().empty())
    {
      // Replaces - We'll sum up the VFs that we can replace in a zone.
      AXOM_PERF_MARK_SECTION("compute_vf_writable", {
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) { vf_writable[i] = 0.; });
        for(const auto& name : shape.getMaterialsReplaced())
        {
          auto mat = getMaterial(name);
          GridFunctionView<ExecSpace> matVFView(mat.first, false);
          axom::for_all<ExecSpace>(
            dataSize,
            AXOM_LAMBDA(axom::IndexType i) {
              vf_writable[i] += matVFView[i];
              vf_writable[i] = (vf_writable[i] > 1.) ? 1. : vf_writable[i];
            });
        }
      });
    }
    else
    {
      // Does not replace. We can replace all except for listed mats.
      AXOM_PERF_MARK_SECTION("compute_vf_writable", {
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) { vf_writable[i] = 1.; });
        for(auto& gf : excludeVFs)
        {
          GridFunctionView<ExecSpace> matVFView(gf, false);
          axom::for_all<ExecSpace>(
            dataSize,
            AXOM_LAMBDA(axom::IndexType i) {
              vf_writable[i] -= matVFView[i];
              vf_writable[i] = (vf_writable[i] < 0.) ? 0. : vf_writable[i];
            });
        }
      });
    }

    // Compute the volume fractions for the current shape's material.
    AXOM_PERF_MARK_SECTION("compute_vf", {
      GridFunctionView<ExecSpace> matVFView(matVF.first);
      GridFunctionView<ExecSpace> shapeVFView(shapeVolFrac);

      // Workaround for HIP so we do not capture through "this" pointer.
      const double* local_overlap_volumes = m_overlap_volumes;
      const double* local_hex_volumes = m_hex_volumes;
      axom::for_all<ExecSpace>(
        dataSize,
        AXOM_LAMBDA(axom::IndexType i) {
          // Update this material's VF and vf_subtract, which is the
          // amount to subtract from the gf's in updateVF.
          double vf = (local_overlap_volumes[i] / local_hex_volumes[i]);

          // Write at most the writable amount.
          double vf_actual = (vf <= vf_writable[i]) ? vf : vf_writable[i];

          // NOTE: if matVFView[i] temporarily exceeds 1, it will be corrected
          //       during the subtraction stage.
          matVFView[i] += vf_actual;
          vf_subtract[i] = vf_actual;

          // Store the max shape VF.
          shapeVFView[i] = vf;
        });
    });

    // Iterate over updateVFs to subtract off VFs we allocated to the
    // current shape's material.
    AXOM_PERF_MARK_SECTION("update_vf", {
      for(auto& gf : updateVFs)
      {
        GridFunctionView<ExecSpace> matVFView(gf);
        axom::for_all<ExecSpace>(
          dataSize,
          AXOM_LAMBDA(axom::IndexType i) {
            constexpr double INSIGNIFICANT_VOLFRAC = 1.e-14;
            double s =
              (matVFView[i] <= vf_subtract[i]) ? matVFView[i] : vf_subtract[i];
            matVFView[i] -= s;
            // Turn any slight negatives or positive insignificant volume fractions to zero.
            matVFView[i] =
              (matVFView[i] < INSIGNIFICANT_VOLFRAC) ? 0. : matVFView[i];
            vf_subtract[i] -= s;
          });
      }
    });
  }
#endif

  /*!
   * \brief Apply material replacement rules for the current shape, using
   *        the appropriate execution policy.
   */
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
  HexahedronType* m_hexes {nullptr};
  BoundingBoxType* m_hex_bbs {nullptr};

  std::vector<mfem::GridFunction*> m_vf_grid_functions;
  std::vector<std::string> m_vf_material_names;
#endif
};

}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_INTERSECTION_SHAPER__HPP_
