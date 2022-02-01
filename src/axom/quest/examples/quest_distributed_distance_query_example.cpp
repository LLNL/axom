// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
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

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"

#include "axom/quest/DistributedClosestPoint.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifndef AXOM_USE_MFEM
  #error This example requires Axom to be configured with MFEM and the AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION option
#endif
#include "mfem.hpp"

#ifndef AXOM_USE_MPI
  #error This example requires Axom to be configured with MPI
#endif
#include "mpi.h"

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
namespace mint = axom::mint;
namespace numerics = axom::numerics;

/// Enum for RAJA runtime policy
enum class RuntimePolicy
{
  seq = 0,
  omp = 1,
  cuda = 2
};

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;

  double circleRadius {1.0};
  int circlePoints {100};
  RuntimePolicy policy {RuntimePolicy::seq};

private:
  bool m_verboseOutput {false};

  // clang-format off
  const std::map<std::string, RuntimePolicy> s_validPolicies{
    #if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
      {"seq", RuntimePolicy::seq}
      #ifdef AXOM_USE_OPENMP
    , {"omp", RuntimePolicy::omp}
      #endif
      #ifdef AXOM_USE_CUDA
    , {"cuda", RuntimePolicy::cuda}
      #endif
    #endif
  };
  // clang-format on

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

    app.add_option("-p, --policy", policy)
      ->description("Set runtime policy for point query method")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(s_validPolicies));

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);
  }
};

/**
 *  \brief Simple wrapper to a blueprint particle mesh
 * 
 *  Given a sidre Group, creates the stubs for a mesh blueptint particle mesh
 */
struct BlueprintParticleMesh
{
public:
  explicit BlueprintParticleMesh(sidre::Group* group = nullptr,
                                 const std::string& coordset = "coords",
                                 const std::string& topology = "mesh")
    : m_group(group)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nranks);

    setBlueprintGroup(m_group, coordset, topology);
  }
  /// Gets the root group for this mesh blueprint
  sidre::Group* rootGroup() const { return m_group; }
  /// Gets the parent group for the blueprint coordinate set
  sidre::Group* coordsGroup() const { return m_coordsGroup; }
  /// Gets the parent group for the blueprint mesh topology
  sidre::Group* topoGroup() const { return m_topoGroup; }

  /// Gets the MPI rank for this mesh
  int getRank() const { return m_rank; }
  /// Gets the number of ranks in the problem
  int getNumRanks() const { return m_nranks; }

  /// Returns true if points have been added to the particle mesh
  bool hasPoints() const
  {
    return m_coordsGroup != nullptr && m_coordsGroup->hasGroup("values");
  }

  /// Returns the number of points in the particle mesh
  int numPoints() const
  {
    return hasPoints() ? m_coordsGroup->getView("values/x")->getNumElements() : 0;
  }

  int dimension() const
  {
    if(hasPoints())
    {
      return m_coordsGroup->hasView("values/z")
        ? 3
        : (m_coordsGroup->hasView("values/y") ? 2 : 1);
    }
    return 0;
  }

  /** 
   * Sets the parent group for the entire mesh and sets up the blueprint stubs
   * for the "coordset", "topologies", "fields" and "state"
   */
  void setBlueprintGroup(sidre::Group* group,
                         const std::string& coordset = "coords",
                         const std::string& topology = "mesh")
  {
    // TODO: Ensure that we delete previous hierarchy if it existed

    m_group = group;

    if(m_group != nullptr)
    {
      createBlueprintStubs(coordset, topology);
    }
  }

  /// Set the coordinate data from an array of primal Points, templated on the dimension
  template <int NDIMS>
  void setPoints(const axom::Array<primal::Point<double, NDIMS>>& pts)
  {
    SLIC_ASSERT_MSG(m_group != nullptr,
                    "Must set blueprint group before setPoints()");

    const int SZ = pts.size();

    // create views into a shared buffer for the coordinates, with stride NDIMS
    auto* buf =
      m_group->getDataStore()->createBuffer(sidre::DOUBLE_ID, NDIMS * SZ)->allocate();
    switch(NDIMS)
    {
    case 3:
      m_coordsGroup->createView("values/z")->attachBuffer(buf)->apply(SZ, 2, NDIMS);
    case 2:  // intentional fallthrough
      m_coordsGroup->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, NDIMS);
    default:  // intentional fallthrough
      m_coordsGroup->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, NDIMS);
    }

    // copy coordinate data into the buffer
    const std::size_t nbytes = sizeof(double) * SZ * NDIMS;
    axom::copy(buf->getVoidPtr(), pts.data(), nbytes);

    // set the default connectivity
    sidre::Array<int> arr(m_topoGroup->createView("elements/connectivity"), SZ, SZ);
    for(int i = 0; i < SZ; ++i)
    {
      arr[i] = i;
    }
  }

  template <typename T>
  void registerNodalScalarField(const std::string& fieldName)
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot register a field with the BlueprintParticleMesh "
                    "before adding points");

    auto* fld = m_fieldsGroup->createGroup(fieldName);
    fld->createViewString("association", "vertex");
    fld->createViewString("topology", m_topoGroup->getName());
    fld->createViewAndAllocate("values",
                               sidre::detail::SidreTT<T>::id,
                               numPoints());
  }

  template <typename T>
  void registerNodalVectorField(const std::string& fieldName)
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot register a field with the BlueprintParticleMesh "
                    "before adding points");

    const int SZ = numPoints();
    const int DIM = dimension();

    auto* fld = m_fieldsGroup->createGroup(fieldName);
    fld->createViewString("association", "vertex");
    fld->createViewString("topology", m_topoGroup->getName());

    // create views into a shared buffer for the coordinates, with stride NDIMS
    auto* buf = m_group->getDataStore()
                  ->createBuffer(sidre::detail::SidreTT<T>::id, DIM * SZ)
                  ->allocate();
    switch(DIM)
    {
    case 3:
      fld->createView("values/z")->attachBuffer(buf)->apply(SZ, 2, DIM);
    case 2:  // intentional fallthrough
      fld->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, DIM);
    default:  // intentional fallthrough
      fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
    }
  }

  bool hasField(const std::string& fieldName) const
  {
    return m_fieldsGroup->hasGroup(fieldName);
  }

  template <typename T>
  axom::ArrayView<T> getNodalScalarField(const std::string& fieldName) const
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot extract a field from the BlueprintParticleMesh "
                    "before adding points");

    T* data = hasField(fieldName)
      ? static_cast<T*>(
          m_fieldsGroup->getView(axom::fmt::format("{}/values", fieldName))
            ->getVoidPtr())
      : nullptr;

    return axom::ArrayView<T>(data, numPoints());
  }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    conduit::Node mesh_node;
    m_group->createNativeLayout(mesh_node);

    bool success = true;
    conduit::Node info;
    if(!conduit::blueprint::mpi::verify("mesh", mesh_node, info, MPI_COMM_WORLD))
    {
      SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
      success = false;
    }

    return success;
  }

  /// Outputs the object mesh to disk
  void saveMesh(const std::string& outputMesh)
  {
    auto* ds = m_group->getDataStore();
    sidre::IOManager writer(MPI_COMM_WORLD);
    writer.write(ds->getRoot(), 1, outputMesh, "sidre_hdf5");

    MPI_Barrier(MPI_COMM_WORLD);

    // Add the bp index to the root file
    writer.writeBlueprintIndexToRootFile(m_group->getDataStore(),
                                         m_group->getPathName(),
                                         outputMesh + ".root",
                                         m_group->getName());
  }

private:
  /// Creates blueprint stubs for this mesh
  void createBlueprintStubs(const std::string& coords, const std::string& topo)
  {
    SLIC_ASSERT(m_group != nullptr);

    m_coordsGroup = m_group->createGroup("coordsets")->createGroup(coords);
    m_coordsGroup->createViewString("type", "explicit");
    m_coordsGroup->createGroup("values");

    m_topoGroup = m_group->createGroup("topologies")->createGroup(topo);
    m_topoGroup->createViewString("coordset", coords);
    m_topoGroup->createViewString("type", "unstructured");
    m_topoGroup->createViewString("elements/shape", "point");

    m_fieldsGroup = m_group->createGroup("fields");

    m_group->createViewScalar<axom::int32>("state/domain_id", m_rank);
  }

private:
  sidre::Group* m_group;

  sidre::Group* m_coordsGroup;
  sidre::Group* m_topoGroup;
  sidre::Group* m_fieldsGroup;

  int m_rank;
  int m_nranks;
};

/** 
 * Helper class to generate a mesh blueprint-conforming particle mesh for the input object.
 * The mesh is represented using a Sidre hierarchy
 */
class ObjectMeshWrapper
{
public:
  ObjectMeshWrapper(sidre::Group* group) : m_group(group), m_mesh(m_group)
  {
    SLIC_ASSERT(m_group != nullptr);
  }

  /// Get a pointer to the root group for this mesh
  sidre::Group* getBlueprintGroup() const { return m_group; }

  std::string getCoordsetName() const
  {
    return m_mesh.coordsGroup()->getName();
  }

  /** 
   * Generates a collection of \a numPoints points along a circle
   * of radius \a radius centered at the origin
   */
  void generateCircleMesh(double radius, int numPoints)
  {
    using axom::utilities::random_real;

    // Check that we're starting with a valid group
    SLIC_ASSERT(m_group != nullptr);

    constexpr int DIM = 2;
    using PointType = primal::Point<double, DIM>;
    using PointArray = axom::Array<PointType>;

    PointArray pts(0, numPoints);
    for(int i = 0; i < numPoints; ++i)
    {
      const double angleInRadians = random_real(0., 2 * M_PI);
      const double rsinT = radius * std::sin(angleInRadians);
      const double rcosT = radius * std::cos(angleInRadians);

      pts.push_back(PointType {rcosT, rsinT});
    }

    m_mesh.setPoints(pts);

    SLIC_ASSERT(m_mesh.isValid());
  }

  /// Outputs the object mesh to disk
  void saveMesh(const std::string& outputMesh = "object_mesh")
  {
    SLIC_INFO(axom::fmt::format(
      "{:=^80}",
      axom::fmt::format("Saving particle mesh '{}' to disk", outputMesh)));

    m_mesh.saveMesh(outputMesh);
  }

private:
  sidre::Group* m_group;
  BlueprintParticleMesh m_mesh;
};

class QueryMeshWrapper
{
public:
  QueryMeshWrapper() : m_dc("closest_point", nullptr, true) { }

  // Returns a pointer to the MFEMSidreDataCollection
  sidre::MFEMSidreDataCollection* getDC() { return &m_dc; }
  const sidre::MFEMSidreDataCollection* getDC() const { return &m_dc; }

  const BlueprintParticleMesh& getParticleMesh() const { return m_queryMesh; }

  sidre::Group* getBlueprintGroup() const { return m_queryMesh.rootGroup(); }

  std::string getCoordsetName() const
  {
    return m_queryMesh.coordsGroup()->getName();
  }

  /// Returns an array containing the positions of the mesh vertices
  template <typename PointArray>
  PointArray getVertexPositions()
  {
    auto* mesh = m_dc.GetMesh();
    const int NV = mesh->GetNV();
    PointArray arr(NV);

    for(auto i : slam::PositionSet<int>(NV))
    {
      mesh->GetNode(i, arr[i].data());
    }

    return arr;
  }

  /// Saves the data collection to disk
  void saveMesh()
  {
    SLIC_INFO(
      axom::fmt::format("{:=^80}",
                        axom::fmt::format("Saving query mesh '{}' to disk",
                                          m_dc.GetCollectionName())));

    m_dc.Save();
  }

  void setupParticleMesh()
  {
    using PointArray2D = axom::Array<primal::Point<double, 2>>;
    using PointArray3D = axom::Array<primal::Point<double, 3>>;

    auto* dsRoot = m_dc.GetBPGroup()->getDataStore()->getRoot();
    m_queryMesh = BlueprintParticleMesh(dsRoot->createGroup("query_mesh"));

    const int DIM = m_dc.GetMesh()->Dimension();

    switch(DIM)
    {
    case 2:
      m_queryMesh.setPoints<2>(getVertexPositions<PointArray2D>());
      break;
    case 3:
      m_queryMesh.setPoints<3>(getVertexPositions<PointArray3D>());
      break;
    }

    m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_rank");
    m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_index");
    m_queryMesh.registerNodalScalarField<double>("min_distance");
    m_queryMesh.registerNodalVectorField<double>("closest_point");

    SLIC_ASSERT(m_queryMesh.isValid());
  }

  /**
   * Loads the mesh as an MFEMSidreDataCollection with 
   * the following fields: "positions", "distances", "directions"
   */
  void setupMesh(const std::string& fileName, const std::string meshFile)
  {
    SLIC_INFO(axom::fmt::format("{:=^80}",
                                axom::fmt::format("Loading '{}' mesh", fileName)));

    sidre::MFEMSidreDataCollection originalMeshDC(fileName, nullptr, false);
    {
      originalMeshDC.SetComm(MPI_COMM_WORLD);
      originalMeshDC.Load(meshFile, "sidre_hdf5");
    }
    SLIC_ASSERT_MSG(originalMeshDC.GetMesh()->Dimension() == 2,
                    "This application currently only supports 2D meshes");
    // TODO: Check order and apply LOR, if necessary

    const int DIM = originalMeshDC.GetMesh()->Dimension();

    // Create the data collection
    mfem::Mesh* cpMesh = nullptr;
    {
      m_dc.SetMeshNodesName("positions");

      auto* pmesh = dynamic_cast<mfem::ParMesh*>(originalMeshDC.GetMesh());
      cpMesh = (pmesh != nullptr) ? new mfem::ParMesh(*pmesh)
                                  : new mfem::Mesh(*originalMeshDC.GetMesh());
      m_dc.SetMesh(cpMesh);
    }

    // Register the distance and direction grid function
    constexpr int order = 1;
    auto* fec = new mfem::H1_FECollection(order, DIM, mfem::BasisType::Positive);
    mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(cpMesh, fec);
    mfem::GridFunction* distances = new mfem::GridFunction(fes);
    distances->MakeOwner(fec);
    m_dc.RegisterField("distance", distances);

    auto* vfec = new mfem::H1_FECollection(order, DIM, mfem::BasisType::Positive);
    mfem::FiniteElementSpace* vfes =
      new mfem::FiniteElementSpace(cpMesh, vfec, DIM);
    mfem::GridFunction* directions = new mfem::GridFunction(vfes);
    directions->MakeOwner(vfec);
    m_dc.RegisterField("direction", directions);
  }

  /// Prints some info about the mesh
  void printMeshInfo()
  {
    switch(m_dc.GetMesh()->Dimension())
    {
    case 2:
      printMeshInfo<2>();
      break;
    case 3:
      printMeshInfo<3>();
      break;
    }
  }

private:
  /**
  * \brief Print some info about the mesh
  *
  * \note In MPI-based configurations, this is a collective call, but only prints on rank 0
  */
  template <int DIM>
  void printMeshInfo()
  {
    mfem::Mesh* mesh = m_dc.GetMesh();

    int myRank = 0;
    int numElements = mesh->GetNE();

    mfem::Vector mins, maxs;

    auto* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);
    if(pmesh != nullptr)
    {
      pmesh->GetBoundingBox(mins, maxs);
      numElements = pmesh->ReduceInt(numElements);
      myRank = pmesh->GetMyRank();
    }
    else
    {
      mesh->GetBoundingBox(mins, maxs);
    }

    if(myRank == 0)
    {
      SLIC_INFO(axom::fmt::format(
        "Mesh has {} elements and (approximate) bounding box {}",
        numElements,
        primal::BoundingBox<double, DIM>(
          primal::Point<double, DIM>(mins.GetData()),
          primal::Point<double, DIM>(maxs.GetData()))));
    }

    slic::flushStreams();
  }

private:
  sidre::MFEMSidreDataCollection m_dc;
  BlueprintParticleMesh m_queryMesh;
};

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  slic::SimpleLogger logger;
  //slic::setAbortOnWarning(true);

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

    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();

    exit(retval);
  }

  constexpr int DIM = 2;

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  using SeqClosestPointQueryType =
    quest::DistributedClosestPoint<2, axom::SEQ_EXEC>;

  #ifdef AXOM_USE_OPENMP
  using OmpClosestPointQueryType =
    quest::DistributedClosestPoint<2, axom::OMP_EXEC>;
  #endif

  #ifdef AXOM_USE_CUDA
  using CudaClosestPointQueryType =
    quest::DistributedClosestPoint<2, axom::CUDA_EXEC<256>>;
  #endif
#endif

  using PointArray = SeqClosestPointQueryType::PointArray;
  using IndexSet = slam::PositionSet<>;

  //---------------------------------------------------------------------------
  // Load/generate object mesh
  //---------------------------------------------------------------------------
  sidre::DataStore objectDS;
  ObjectMeshWrapper object_mesh_wrapper(
    objectDS.getRoot()->createGroup("object_mesh"));

  object_mesh_wrapper.generateCircleMesh(params.circleRadius,
                                         params.circlePoints);
  object_mesh_wrapper.saveMesh();

  //---------------------------------------------------------------------------
  // Load computational mesh and generate a particle mesh over its nodes
  // These will be used to query the closest points on the object mesh(es)
  //---------------------------------------------------------------------------
  QueryMeshWrapper query_mesh_wrapper;

  query_mesh_wrapper.setupMesh(params.getDCMeshName(), params.meshFile);
  query_mesh_wrapper.printMeshInfo();
  query_mesh_wrapper.setupParticleMesh();

  // Copy mesh nodes into qpts array
  auto qPts = query_mesh_wrapper.getVertexPositions<PointArray>();
  const int nQueryPts = qPts.size();

  // Create an array to hold the object points
  PointArray objectPts;

  //---------------------------------------------------------------------------
  // Initialize spatial index for querying points, and run query
  //---------------------------------------------------------------------------

  auto init_str =
    axom::fmt::format("{:=^80}",
                      axom::fmt::format("Initializing BVH tree over {} points",
                                        params.circlePoints));

  auto query_str = axom::fmt::format(
    "{:=^80}",
    axom::fmt::format("Computing closest points for {} query points", nQueryPts));

  axom::utilities::Timer initTimer(false);
  axom::utilities::Timer queryTimer(false);

  switch(params.policy)
  {
  case RuntimePolicy::seq:
  {
    SeqClosestPointQueryType query;
    query.setVerbosity(params.isVerbose());
    SLIC_INFO(init_str);
    query.setObjectMesh(object_mesh_wrapper.getBlueprintGroup(),
                        object_mesh_wrapper.getCoordsetName());

    SLIC_INFO(init_str);
    initTimer.start();
    query.generateBVHTree();
    initTimer.stop();

    SLIC_INFO(query_str);
    queryTimer.start();
    query.computeClosestPoints(query_mesh_wrapper.getBlueprintGroup(),
                               query_mesh_wrapper.getCoordsetName());
    queryTimer.stop();
    objectPts = query.points();
  }
  break;
  case RuntimePolicy::omp:
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(AXOM_USE_OPENMP)
  {
    OmpClosestPointQueryType query;
    query.setVerbosity(params.isVerbose());
    query.setObjectMesh(object_mesh_wrapper.getBlueprintGroup(),
                        object_mesh_wrapper.getCoordsetName());

    SLIC_INFO(init_str);
    initTimer.start();
    query.generateBVHTree();
    initTimer.stop();

    SLIC_INFO(query_str);
    queryTimer.start();
    query.computeClosestPoints(query_mesh_wrapper.getBlueprintGroup(),
                               query_mesh_wrapper.getCoordsetName());
    queryTimer.stop();

    objectPts = query.points();
  }
#endif
  break;
  case RuntimePolicy::cuda:
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA)
  {
    CudaClosestPointQueryType query;
    query.setVerbosity(params.isVerbose());
    query.setObjectMesh(object_mesh_wrapper.getBlueprintGroup(),
                        object_mesh_wrapper.getCoordsetName());

    SLIC_INFO(init_str);
    initTimer.start();
    query.generateBVHTree();
    initTimer.stop();

    SLIC_INFO(query_str);
    queryTimer.start();
    query.computeClosestPoints(query_mesh_wrapper.getBlueprintGroup(),
                               query_mesh_wrapper.getCoordsetName());
    queryTimer.stop();

    objectPts = query.points();
  }
#endif
  break;
  }

  auto cpIndices =
    query_mesh_wrapper.getParticleMesh().getNodalScalarField<axom::IndexType>(
      "cp_index");

  if(params.isVerbose())
  {
    SLIC_INFO(axom::fmt::format("Closest points ({}):", cpIndices.size()));
    for(auto i : IndexSet(cpIndices.size()))
    {
      SLIC_INFO(axom::fmt::format("\t{}: {}", i, cpIndices[i]));
    }
  }

  SLIC_INFO(axom::fmt::format("Initialization with policy {} took {} seconds",
                              params.policy,
                              initTimer.elapsedTimeInSec()));
  SLIC_INFO(axom::fmt::format("Query with policy {} took {} seconds",
                              params.policy,
                              queryTimer.elapsedTimeInSec()));

  //---------------------------------------------------------------------------
  // Transform closest points to distances and directions
  //---------------------------------------------------------------------------
  using primal::squared_distance;

  auto* distances = query_mesh_wrapper.getDC()->GetField("distance");
  auto* directions = query_mesh_wrapper.getDC()->GetField("direction");
  SLIC_INFO(axom::fmt::format(" distance size: {}", distances->Size()));
  mfem::Array<int> dofs;
  for(auto i : IndexSet(nQueryPts))
  {
    const auto& cp = objectPts[cpIndices[i]];
    (*distances)(i) = sqrt(squared_distance(qPts[i], cp));

    primal::Vector<double, DIM> dir(qPts[i], cp);
    directions->FESpace()->GetVertexVDofs(i, dofs);
    directions->SetSubVector(dofs, dir.data());
  }

  //---------------------------------------------------------------------------
  // Cleanup, save mesh/fields and exit
  //---------------------------------------------------------------------------
  query_mesh_wrapper.saveMesh();

  MPI_Finalize();

  return 0;
}
