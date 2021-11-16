// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_distributed_distance_query_example.cpp
 * \brief Driver for a distributed distance query
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/quest.hpp"
#include "axom/slam.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifndef AXOM_USE_MFEM
  #error This example requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif
#include "mfem.hpp"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// RAJA
#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

// C/C++ includes
#include <string>
#include <map>
#include <cmath>

namespace quest = axom::quest;
namespace slic = axom::slic;
namespace sidre = axom::sidre;
namespace slam = axom::slam;
namespace spin = axom::spin;
namespace primal = axom::primal;
namespace numerics = axom::numerics;

template <int NDIMS = 2, typename ExecSpace = axom::SEQ_EXEC>
class ClosestPointQuery
{
public:
  // TODO: generalize to 3D
  static_assert(NDIMS == 2, "ClosestPointQuery only currently supports 2D");

  static constexpr int DIM = NDIMS;
  using PointType = primal::Point<double, DIM>;
  using BoxType = primal::BoundingBox<double, DIM>;
  using PointArray = axom::Array<PointType>;
  using BoxArray = axom::Array<BoxType>;
  using BVHTreeType = spin::BVH<DIM, ExecSpace>;

private:
  struct MinCandidate
  {
    /// Squared distance to query point
    double minSqDist {numerics::floating_point_limits<double>::max()};
    /// Index within mesh of closest element
    int minElem;
  };

public:
  ClosestPointQuery(
    int allocatorID = axom::execution_space<ExecSpace>::allocatorID())
    : m_allocatorID(allocatorID)
  { }

  /// Utility function to generate an array of 2D points
  void generatePoints(double radius, int numPoints)
  {
    using axom::utilities::random_real;

    m_points.clear();
    m_points.reserve(numPoints);

    for(int i = 0; i < numPoints; ++i)
    {
      const double angleInRadians = random_real(0., 2 * M_PI);
      const double rsinT = radius * std::sin(angleInRadians);
      const double rcosT = radius * std::cos(angleInRadians);

      m_points.emplace_back(Point2D {rcosT, rsinT});
    }
  }

  bool generateBVHTree()
  {
    const int npts = m_points.size();
    BoxType* boxes = axom::allocate<BoxType>(npts, m_allocatorID);
    axom::for_all<ExecSpace>(
      npts,
      AXOM_LAMBDA(axom::IndexType i) { boxes[i] = BoxType {m_points[i]}; });

    // Build bounding volume hierarchy
    m_bvh.setAllocatorID(m_allocatorID);
    int result = m_bvh.initialize(boxes, npts);

    axom::deallocate(boxes);
    return (result == spin::BVH_BUILD_OK);
  }

  void computeClosestPoints(const PointArray& queryPts,
                            axom::Array<axom::IndexType>& cpIndexes) const
  {
    SLIC_ASSERT(!queryPts.empty());

    const int nPts = queryPts.size();

    cpIndexes.resize(nPts);
    axom::IndexType* indexData = cpIndexes.data();

    // Get a device-useable iterator
    auto it = m_bvh.getTraverser();

    using axom::primal::squared_distance;
    using int32 = axom::int32;

    AXOM_PERF_MARK_SECTION(
      "ComputeClosestPoints",
      axom::for_all<ExecSpace>(
        nPts,
        AXOM_LAMBDA(int32 idx) {
          PointType qpt = queryPts[idx];

          MinCandidate curr_min {};

          auto searchMinDist = [&](int32 current_node, const int32* leaf_nodes) {
            const int candidate_idx = leaf_nodes[current_node];
            const PointType candidate_pt = m_points[candidate_idx];
            const double sq_dist = squared_distance(qpt, candidate_pt);

            if(sq_dist < curr_min.minSqDist)
            {
              curr_min.minSqDist = sq_dist;
              curr_min.minElem = candidate_idx;
            }
          };

          auto traversePredicate = [&](const PointType& p,
                                       const BoxType& bb) -> bool {
            return squared_distance(p, bb) <= curr_min.minSqDist;
          };

          // Traverse the tree, searching for the point with minimum distance.
          it.traverse_tree(qpt, searchMinDist, traversePredicate);

          indexData[idx] = curr_min.minElem;
        }););
  }

  const PointArray& points() const { return m_points; }

private:
  PointArray m_points;
  BoxArray m_boxes;
  BVHTreeType m_bvh;

  int m_allocatorID;
};

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;

  double circleRadius {1.0};
  int circlePoints {100};

private:
  bool m_verboseOutput {false};

public:
  bool isVerbose() const { return m_verboseOutput; }

  std::string getDCMeshName() const
  {
    using axom::utilities::string::removeSuffix;

    // Remove the parent directories and file suffix
    std::string name = axom::Path(meshFile).baseName();
    name = removeSuffix(name, ".root");

    return name;
  }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-m,--mesh-file", meshFile)
      ->description(
        "Path to computational mesh (generated by MFEMSidreDataCollection)")
      ->check(axom::CLI::ExistingFile)
      ->required();

    app.add_flag("-v,--verbose,!--no-verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("-r,--radius", circleRadius)
      ->description("Radius for circle")
      ->capture_default_str();

    app.add_option("-n,--num-samples", circlePoints)
      ->description("Number of points for circle")
      ->capture_default_str();

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);
  }
};

/**
 * \brief Print some info about the mesh
 *
 * \note In MPI-based configurations, this is a collective call, but only prints on rank 0
 */
void printMeshInfo(mfem::Mesh* mesh, const std::string& prefixMessage = "")
{
  namespace primal = axom::primal;

  int myRank = 0;
#ifdef AXOM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

  int numElements = mesh->GetNE();

  mfem::Vector mins, maxs;
#ifdef MFEM_USE_MPI
  auto* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);
  if(pmesh != nullptr)
  {
    pmesh->GetBoundingBox(mins, maxs);
    numElements = pmesh->ReduceInt(numElements);
    myRank = pmesh->GetMyRank();
  }
  else
#endif
  {
    mesh->GetBoundingBox(mins, maxs);
  }

  if(myRank == 0)
  {
    switch(mesh->Dimension())
    {
    case 2:
      SLIC_INFO(axom::fmt::format(
        "{} mesh has {} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 2>(primal::Point<double, 2>(mins.GetData()),
                                       primal::Point<double, 2>(maxs.GetData()))));
      break;
    case 3:
      SLIC_INFO(axom::fmt::format(
        "{} mesh has {} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 3>(primal::Point<double, 3>(mins.GetData()),
                                       primal::Point<double, 3>(maxs.GetData()))));
      break;
    }
  }

  slic::flushStreams();
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#else
  int my_rank = 0;
  int num_ranks = 1;
#endif

  slic::SimpleLogger logger;

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  Input params;
  axom::CLI::App app {"Driver for distributed distance query"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(my_rank == 0)
    {
      retval = app.exit(e);
    }

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(retval);
  }

  //---------------------------------------------------------------------------
  // Load mesh
  //---------------------------------------------------------------------------
  constexpr int DIM = 2;

  SLIC_INFO(axom::fmt::format(
    "{:=^80}",
    axom::fmt::format("Loading '{}' mesh", params.getDCMeshName())));

  const bool dc_owns_data = true;
  mfem::Mesh* originalMesh = nullptr;
  sidre::MFEMSidreDataCollection originalMeshDC(params.getDCMeshName(),
                                                originalMesh,
                                                dc_owns_data);
  {
    originalMeshDC.SetComm(MPI_COMM_WORLD);
    std::string protocol = "sidre_hdf5";
    originalMeshDC.Load(params.meshFile, protocol);
  }
  SLIC_ASSERT_MSG(originalMeshDC.GetMesh()->Dimension() == DIM,
                  "This application currently only supports 2D meshes");
  // TODO: Check order and apply LOR, if necessary

  mfem::Mesh* cpMesh = nullptr;
  sidre::MFEMSidreDataCollection closestPointDC("closest_point",
                                                cpMesh,
                                                dc_owns_data);
  {
    closestPointDC.SetMeshNodesName("positions");

    auto* pmesh = dynamic_cast<mfem::ParMesh*>(originalMeshDC.GetMesh());
    cpMesh = (pmesh != nullptr) ? new mfem::ParMesh(*pmesh)
                                : new mfem::Mesh(*originalMeshDC.GetMesh());
    closestPointDC.SetMesh(cpMesh);
  }
  printMeshInfo(closestPointDC.GetMesh(), "After loading");

  //---------------------------------------------------------------------------
  // Initialize spatial index
  //---------------------------------------------------------------------------

  SLIC_INFO(
    axom::fmt::format("{:=^80}",
                      axom::fmt::format("Initializing BVH tree over {} points",
                                        params.circlePoints)));

  using ExecSpace = axom::SEQ_EXEC;
  using ClosestPointQueryType = ClosestPointQuery<2, ExecSpace>;
  using PointArray = ClosestPointQueryType::PointArray;
  ClosestPointQueryType query;

  query.generatePoints(params.circleRadius, params.circlePoints);

  if(params.isVerbose())
  {
    SLIC_INFO("Points on object:");
    const auto& arr = query.points();
    for(auto i : slam::PositionSet<>(arr.size()))
    {
      const double mag = sqrt(arr[i][0] * arr[i][0] + arr[i][1] * arr[i][1]);
      SLIC_INFO(axom::fmt::format("\t{}: {} -- {}", i, arr[i], mag));
    }
  }

  // Generate BVH tree over the points
  query.generateBVHTree();

  //---------------------------------------------------------------------------
  // Set up query points
  //---------------------------------------------------------------------------
  const int nQueryPts = cpMesh->GetNV();

  SLIC_INFO(axom::fmt::format(
    "{:=^80}",
    axom::fmt::format("Computing closest points for {} query points", nQueryPts)));

  PointArray qPts(nQueryPts);
  for(auto i : slam::PositionSet<>(nQueryPts))
  {
    cpMesh->GetNode(i, qPts[i].data());
  }

  // Register the distances grid function
  constexpr int order = 1;
  auto* fec = new mfem::H1_FECollection(order, DIM, mfem::BasisType::Positive);
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(cpMesh, fec);
  mfem::GridFunction* distances = new mfem::GridFunction(fes);
  distances->MakeOwner(fec);
  closestPointDC.RegisterField("distance", distances);

  using IndexArray = axom::Array<axom::IndexType>;
  IndexArray cpIndices;
  query.computeClosestPoints(qPts, cpIndices);

  if(params.isVerbose())
  {
    SLIC_INFO(axom::fmt::format("Closest points ({}):", cpIndices.size()));
    for(auto i : slam::PositionSet<>(cpIndices.size()))
    {
      SLIC_INFO(axom::fmt::format("\t{}: {}", i, cpIndices[i]));
    }
  }

  SLIC_INFO(axom::fmt::format(" distance size: {}", distances->Size()));

  using primal::squared_distance;
  const auto& objectPts = query.points();
  for(auto i : slam::PositionSet<>(nQueryPts))
  {
    (*distances)[i] = sqrt(squared_distance(qPts[i], objectPts[cpIndices[i]]));
  }

  //---------------------------------------------------------------------------
  // Save meshes and fields
  //---------------------------------------------------------------------------
#ifdef MFEM_USE_MPI
  SLIC_INFO(
    axom::fmt::format("{:=^80}",
                      axom::fmt::format("Saving mesh '{}' to disk",
                                        closestPointDC.GetCollectionName())));

  closestPointDC.Save();
#endif

  //---------------------------------------------------------------------------
  // Cleanup and exit
  //---------------------------------------------------------------------------

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
