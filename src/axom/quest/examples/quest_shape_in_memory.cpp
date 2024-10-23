// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_shape_mint_mesh.cpp
 * \brief Driver application for shaping material volume fractions onto a simulation mesh
 * using a mint::UnstructuredMesh of tets.
 * Modeled after shaping_driver.cpp.
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/klee.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#if !defined(AXOM_USE_CONDUIT)
  #error Shaping functionality requires Axom to be configured with Conduit
#endif

#include "conduit_relay_io_blueprint.hpp"

#if defined(AXOM_USE_MFEM)
  #include "mfem.hpp"
#endif

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

// RAJA
#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

// C/C++ includes
#include <string>
#include <vector>
#include <memory>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;
namespace sidre = axom::sidre;

//------------------------------------------------------------------------------

using RuntimePolicy = axom::runtime_policy::Policy;

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string outputFile;

  // Values for some shape geometries
  // See createShape_*() for how specific shapes use them.
  double radius {1.0};
  double radius2 {0.3};
  double length {2.0};
  std::vector<double> center {0.1, 0.2, 0.3};
  std::vector<double> direction {0.1, 0.2, 0.4};

  // Shape transformation parameters
  std::vector<double> scaleFactors;

  // Mesh format: mfem or blueprint
  std::string meshFormat = "mfem";

  // Inline mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int boxDim {-1};

  // The shape to run.
  std::string testShape {"tetmesh"};
  // The shapes this example is set up to run.
  const std::set<std::string> availableShapes {"tetmesh",
                                               "sphere",
                                               "cyl",
                                               "cone",
                                               "vor",
                                               "tet",
                                               "hex",
                                               "plane"};

  RuntimePolicy policy {RuntimePolicy::seq};
  int quadratureOrder {5};
  int outputOrder {2};
  int refinementLevel {7};
  double weldThresh {1e-9};
  double percentError {-1.};
  std::string annotationMode {"none"};

  std::string backgroundMaterial;

  // clang-format off
  enum class MeshType { bp = 0, mfem = 1 };
  const std::map<std::string, MeshType> meshTypeChoices
    { {"bp", MeshType::bp} , {"mfem", MeshType::mfem} };
  // clang-format on
  MeshType meshType {MeshType::bp};
  bool useMfem() { return meshType == MeshType::mfem; }
  bool useBlueprint() { return meshType == MeshType::bp; }

private:
  bool m_verboseOutput {false};

public:
  bool isVerbose() const { return m_verboseOutput; }

  /// @brief Return volume of input box mesh
  double boxMeshVolume() const
  {
    primal::Vector<double, 3> x {boxMaxs[0] - boxMins[0], 0, 0};
    primal::Vector<double, 3> y {0, boxMaxs[1] - boxMins[1], 0};
    primal::Vector<double, 3> z {0, 0, boxMaxs[2] - boxMins[2]};
    double volume = primal::Vector<double, 3>::scalar_triple_product(x, y, z);
    return volume;
  }

#if defined(AXOM_USE_MFEM)
  /// Generate an mfem Cartesian mesh, scaled to the bounding box range
  mfem::Mesh* createBoxMesh()
  {
    mfem::Mesh* mesh = nullptr;

    switch(boxDim)
    {
    case 2:
    {
      using BBox2D = primal::BoundingBox<double, 2>;
      using Pt2D = primal::Point<double, 2>;
      auto res = primal::NumericArray<int, 2>(boxResolution.data());
      auto bbox = BBox2D(Pt2D(boxMins.data()), Pt2D(boxMaxs.data()));

      SLIC_INFO(axom::fmt::format(
        "Creating inline box mesh of resolution {} and bounding box {}",
        res,
        bbox));

      mesh = quest::util::make_cartesian_mfem_mesh_2D(bbox, res, outputOrder);
    }
    break;
    case 3:
    {
      using BBox3D = primal::BoundingBox<double, 3>;
      using Pt3D = primal::Point<double, 3>;
      auto res = primal::NumericArray<int, 3>(boxResolution.data());
      auto bbox = BBox3D(Pt3D(boxMins.data()), Pt3D(boxMaxs.data()));

      SLIC_INFO(axom::fmt::format(
        "Creating inline box mesh of resolution {} and bounding box {}",
        res,
        bbox));

      mesh = quest::util::make_cartesian_mfem_mesh_3D(bbox, res, outputOrder);
    }
    break;
    default:
      SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
      break;
    }

    // Handle conversion to parallel mfem mesh
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    {
      int* partitioning = nullptr;
      int part_method = 0;
      mfem::Mesh* parallelMesh =
        new mfem::ParMesh(MPI_COMM_WORLD, *mesh, partitioning, part_method);
      delete[] partitioning;
      delete mesh;
      mesh = parallelMesh;
    }
  #endif

    return mesh;
  }

  std::unique_ptr<sidre::MFEMSidreDataCollection> loadComputationalMesh()
  {
    constexpr bool dc_owns_data = true;
    mfem::Mesh* mesh = createBoxMesh();
    std::string name = "mesh";

    auto dc = std::unique_ptr<sidre::MFEMSidreDataCollection>(
      new sidre::MFEMSidreDataCollection(name, mesh, dc_owns_data));
    dc->SetComm(MPI_COMM_WORLD);

    return dc;
  }
#endif

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-o,--outputFile", outputFile)
      ->description("Path to output file(s)");

    app.add_flag("-v,--verbose,!--no-verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("--meshType", meshType)
      ->description("Type of computational mesh to shape on")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(meshTypeChoices));

    app.add_option("-t,--weld-threshold", weldThresh)
      ->description("Threshold for welding")
      ->check(axom::CLI::NonNegativeNumber)
      ->capture_default_str();

    app.add_option("-e,--percent-error", percentError)
      ->description(
        "Percent error used for calculating curve refinement and revolved "
        "volume.\n"
        "If this value is provided then dynamic curve refinement will be used\n"
        "instead of segment-based curve refinement.")
      ->check(axom::CLI::PositiveNumber)
      ->capture_default_str();

    app.add_option("-s,--testShape", testShape)
      ->description("The shape to run")
      ->check(axom::CLI::IsMember(availableShapes));

#ifdef AXOM_USE_CALIPER
    app.add_option("--caliper", annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        "Use 'help' to see full list.")
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    app.add_option("--center", center)
      ->description("Center of sphere or base of cone/cyl/VOR (x,y[,z]) shape")
      ->expected(2, 3);

    app.add_option("--radius", radius)
      ->description("Radius of sphere or cylinder shape")
      ->check(axom::CLI::PositiveNumber);

    app.add_option("--length", length)
      ->description("Length of cone/cyl/VOR shape")
      ->check(axom::CLI::PositiveNumber);

    app.add_option("--dir", direction)
      ->description("Direction of axis of cone/cyl/VOR (x,y[,z]) shape")
      ->expected(2, 3);

    app.add_option("--radius2", radius2)
      ->description("Second radius of cone shape")
      ->check(axom::CLI::PositiveNumber);

    app.add_option("--scale", scaleFactors)
      ->description("Scale factor to apply to shape (x,y[,z])")
      ->expected(2, 3)
      ->check(axom::CLI::PositiveNumber);

    // use either an input mesh file or a simple inline Cartesian mesh
    {
      auto* inline_mesh_subcommand =
        app.add_subcommand("inline_mesh")
          ->description("Options for setting up a simple inline mesh")
          ->fallthrough();

      inline_mesh_subcommand->add_option("--min", boxMins)
        ->description("Min bounds for box mesh (x,y[,z])")
        ->expected(2, 3)
        ->required();
      inline_mesh_subcommand->add_option("--max", boxMaxs)
        ->description("Max bounds for box mesh (x,y[,z])")
        ->expected(2, 3)
        ->required();

      inline_mesh_subcommand->add_option("--res", boxResolution)
        ->description("Resolution of the box mesh (i,j[,k])")
        ->expected(2, 3)
        ->required();

      inline_mesh_subcommand->add_option("-d,--dimension", boxDim)
        ->description("Dimension of the box mesh")
        ->check(axom::CLI::PositiveNumber)
        ->required();
    }

    app.add_option("--background-material", backgroundMaterial)
      ->description("Sets the name of the background material");

    // parameters that only apply to the intersection method
    {
      auto* intersection_options =
        app.add_option_group("intersection",
                             "Options related to intersection-based queries");

      intersection_options->add_option("-r, --refinements", refinementLevel)
        ->description("Number of refinements to perform for revolved contour")
        ->capture_default_str()
        ->check(axom::CLI::NonNegativeNumber);

      std::stringstream pol_sstr;
      pol_sstr << "Set runtime policy for intersection-based sampling method.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
      pol_sstr << "\nSet to 'seq' or 0 to use the RAJA sequential policy.";
  #ifdef AXOM_USE_OPENMP
      pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
  #endif
  #ifdef AXOM_USE_CUDA
      pol_sstr << "\nSet to 'cuda' or 2 to use the RAJA CUDA policy.";
  #endif
  #ifdef AXOM_USE_HIP
      pol_sstr << "\nSet to 'hip' or 3 to use the RAJA HIP policy.";
  #endif
#endif

      intersection_options->add_option("-p, --policy", policy, pol_sstr.str())
        ->capture_default_str()
        ->transform(
          axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));
    }
    app.get_formatter()->column_width(50);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);
  }
};  // struct Input
Input params;

#if defined(AXOM_USE_MFEM)
/**
 * \brief Print some info about the mesh
 *
 * \note In MPI-based configurations, this is a collective call, but
 * only prints on rank 0
 */
void printMfemMeshInfo(mfem::Mesh* mesh, const std::string& prefixMessage = "")
{
  namespace primal = axom::primal;

  int myRank = 0;
  #ifdef AXOM_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  #endif

  int numElements = mesh->GetNE();

  mfem::Vector mins, maxs;
  #ifdef MFEM_USE_MPI
  auto* parallelMesh = dynamic_cast<mfem::ParMesh*>(mesh);
  if(parallelMesh != nullptr)
  {
    parallelMesh->GetBoundingBox(mins, maxs);
    numElements = parallelMesh->ReduceInt(numElements);
    myRank = parallelMesh->GetMyRank();
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
        axom::utilities::locale(),
        "{} mesh has {:L} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 2>(primal::Point<double, 2>(mins.GetData()),
                                       primal::Point<double, 2>(maxs.GetData()))));
      break;
    case 3:
      SLIC_INFO(axom::fmt::format(
        axom::utilities::locale(),
        "{} mesh has {:L} elements and (approximate) bounding box {}",
        prefixMessage,
        numElements,
        primal::BoundingBox<double, 3>(primal::Point<double, 3>(mins.GetData()),
                                       primal::Point<double, 3>(maxs.GetData()))));
      break;
    }
  }

  slic::flushStreams();
}
#endif

const std::string topoName = "mesh";
const std::string coordsetName = "coords";

axom::sidre::Group* createBoxMesh(axom::sidre::Group* meshGrp)
{
  using BBox3D = primal::BoundingBox<double, 3>;
  using Pt3D = primal::Point<double, 3>;
  auto res = primal::NumericArray<int, 3>(params.boxResolution.data());
  auto bbox = BBox3D(Pt3D(params.boxMins.data()), Pt3D(params.boxMaxs.data()));
  axom::quest::util::make_unstructured_blueprint_box_mesh(meshGrp,
                                                          bbox,
                                                          res,
                                                          topoName,
                                                          coordsetName);
#if defined(AXOM_DEBUG)
  conduit::Node meshNode, info;
  meshGrp->createNativeLayout(meshNode);
  SLIC_ASSERT(conduit::blueprint::mesh::verify(meshNode, info));
#endif

  return meshGrp;
}

/// \brief Utility function to initialize the logger
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream {nullptr};

#ifdef AXOM_USE_MPI
  int num_ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  if(num_ranks > 1)
  {
    std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
    const int RLIMIT = 8;
    logStream =
      new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
    logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
  }
  else
#endif  // AXOM_USE_MPI
  {
    std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
    logStream = new slic::GenericOutputStream(&std::cout, fmt);
  }

  slic::addStreamToAllMsgLevels(logStream);
}

/// \brief Utility function to finalize the logger
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

// Single triangle ShapeSet.
axom::klee::ShapeSet create2DShapeSet(sidre::DataStore& ds)
{
  sidre::Group* meshGroup = ds.getRoot()->createGroup("triangleMesh");
  AXOM_UNUSED_VAR(meshGroup);  // variable is only referenced in debug configs
  const std::string topo = "mesh";
  const std::string coordset = "coords";
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> triangleMesh(
    2,
    axom::mint::CellType::TRIANGLE,
    meshGroup,
    topo,
    coordset);

  double lll = 2.0;

  // Insert tet at origin.
  triangleMesh.appendNode(0.0, 0.0);
  triangleMesh.appendNode(lll, 0.0);
  triangleMesh.appendNode(0.0, lll);
  axom::IndexType conn[3] = {0, 1, 2};
  triangleMesh.appendCell(conn);

  SLIC_ASSERT(axom::mint::blueprint::isValidRootGroup(meshGroup));

  axom::klee::TransformableGeometryProperties prop {
    axom::klee::Dimensions::Two,
    axom::klee::LengthUnit::unspecified};
  axom::klee::Geometry triangleGeom(prop,
                                    triangleMesh.getSidreGroup(),
                                    topo,
                                    nullptr);

  std::vector<axom::klee::Shape> shapes;
  axom::klee::Shape triangleShape(
    "triangle",
    "AL",
    {},
    {},
    axom::klee::Geometry {prop, triangleMesh.getSidreGroup(), topo, nullptr});
  shapes.push_back(axom::klee::Shape {
    "triangle",
    "AL",
    {},
    {},
    axom::klee::Geometry {prop, triangleMesh.getSidreGroup(), topo, nullptr}});

  axom::klee::ShapeSet shapeSet;
  shapeSet.setShapes(shapes);
  shapeSet.setDimensions(axom::klee::Dimensions::Two);

  return shapeSet;
}

axom::klee::Shape createShape_Sphere()
{
  axom::primal::Sphere<double, 3> sphere {params.center.data(), params.radius};

  axom::klee::TransformableGeometryProperties prop {
    axom::klee::Dimensions::Three,
    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  std::shared_ptr<axom::klee::Scale> scaleOp;
  if(!params.scaleFactors.empty())
  {
    scaleOp = std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                                  params.scaleFactors[1],
                                                  params.scaleFactors[2],
                                                  prop);
  }

  const axom::IndexType levelOfRefinement = params.refinementLevel;
  axom::klee::Geometry sphereGeometry(prop, sphere, levelOfRefinement, scaleOp);
  axom::klee::Shape sphereShape("sphere", "AU", {}, {}, sphereGeometry);

  return sphereShape;
}

axom::klee::Shape createShape_TetMesh(sidre::DataStore& ds)
{
  // Shape a single tetrahedron.
  sidre::Group* meshGroup = ds.getRoot()->createGroup("tetMesh");
  AXOM_UNUSED_VAR(meshGroup);  // variable is only referenced in debug configs
  const std::string topo = "mesh";
  const std::string coordset = "coords";
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> tetMesh(
    3,
    axom::mint::CellType::TET,
    meshGroup,
    topo,
    coordset);

  double lll = params.length;
  double h = lll / 2;  // To center the shape.

  tetMesh.appendNode(0.0 - h, 0.0 - h, 0.0 - h);
  tetMesh.appendNode(lll - h, 0.0 - h, 0.0 - h);
  tetMesh.appendNode(0.0 - h, lll - h, 0.0 - h);
  tetMesh.appendNode(0.0 - h, 0.0 - h, lll - h);
  tetMesh.appendNode(lll - h, lll - h, lll - h);
  tetMesh.appendNode(lll - h, lll - h, 0.0 - h);
  tetMesh.appendNode(0.0 - h, lll - h, lll - h);
  tetMesh.appendNode(lll - h, 0.0 - h, lll - h);
  axom::IndexType conn0[4] = {0, 1, 2, 3};
  tetMesh.appendCell(conn0);
  axom::IndexType conn1[4] = {4, 5, 6, 7};
  tetMesh.appendCell(conn1);

  SLIC_ASSERT(axom::mint::blueprint::isValidRootGroup(meshGroup));

  axom::klee::TransformableGeometryProperties prop {
    axom::klee::Dimensions::Three,
    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  std::shared_ptr<axom::klee::Scale> scaleOp;
  if(!params.scaleFactors.empty())
  {
    scaleOp = std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                                  params.scaleFactors[1],
                                                  params.scaleFactors[2],
                                                  prop);
  }

  axom::klee::Geometry tetMeshGeometry(prop,
                                       tetMesh.getSidreGroup(),
                                       topo,
                                       {scaleOp});
  axom::klee::Shape tetShape("tetmesh", "TETMESH", {}, {}, tetMeshGeometry);

  return tetShape;
}

axom::klee::Geometry createGeometry_Vor(axom::primal::Point<double, 3>& vorBase,
                                        axom::primal::Vector<double, 3>& vorDirection,
                                        axom::Array<double, 2>& discreteFunction)
{
  axom::klee::TransformableGeometryProperties prop {
    axom::klee::Dimensions::Three,
    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  std::shared_ptr<axom::klee::Scale> scaleOp;
  if(!params.scaleFactors.empty())
  {
    scaleOp = std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                                  params.scaleFactors[1],
                                                  params.scaleFactors[2],
                                                  prop);
  }

  const axom::IndexType levelOfRefinement = params.refinementLevel;
  axom::klee::Geometry vorGeometry(prop,
                                   discreteFunction,
                                   vorBase,
                                   vorDirection,
                                   levelOfRefinement,
                                   scaleOp);
  return vorGeometry;
}

axom::klee::Shape createShape_Vor()
{
  Point3D vorBase {0.0, 0.0, -params.length / 2};
  axom::primal::Vector<double, 3> vorDirection {params.direction.data()};
  int numIntervals = 5;
  axom::Array<double, 2> discreteFunction({numIntervals + 1, 2},
                                          axom::ArrayStrideOrder::ROW);
  double dz = params.length / numIntervals;
  discreteFunction[0][0] = 0 * dz;
  discreteFunction[0][1] = params.radius;
  discreteFunction[1][0] = 1 * dz;
  discreteFunction[1][1] = params.radius;
  discreteFunction[2][0] = 2 * dz;
  discreteFunction[2][1] = params.radius2;
  discreteFunction[3][0] = 3 * dz;
  discreteFunction[3][1] = params.radius2;
  discreteFunction[4][0] = 4 * dz;
  discreteFunction[4][1] = params.radius;
  discreteFunction[5][0] = 5 * dz;
  discreteFunction[5][1] = 0.0;

  axom::klee::Geometry vorGeometry =
    createGeometry_Vor(vorBase, vorDirection, discreteFunction);

  axom::klee::Shape vorShape("vor", "VOR", {}, {}, vorGeometry);

  return vorShape;
}

axom::klee::Shape createShape_Cylinder()
{
  Point3D vorBase {0.0, 0.0, -params.length / 2};
  axom::primal::Vector<double, 3> vorDirection {params.direction.data()};
  axom::Array<double, 2> discreteFunction({2, 2}, axom::ArrayStrideOrder::ROW);
  double radius = params.radius;
  double height = params.length;
  discreteFunction[0][0] = 0.0;
  discreteFunction[0][1] = radius;
  discreteFunction[1][0] = height;
  discreteFunction[1][1] = radius;

  axom::klee::Geometry vorGeometry =
    createGeometry_Vor(vorBase, vorDirection, discreteFunction);

  axom::klee::Shape vorShape("cyl", "CYL", {}, {}, vorGeometry);

  return vorShape;
}

axom::klee::Shape createShape_Cone()
{
  Point3D vorBase {0.0, 0.0, -params.length / 2};
  axom::primal::Vector<double, 3> vorDirection {params.direction.data()};
  axom::Array<double, 2> discreteFunction({2, 2}, axom::ArrayStrideOrder::ROW);
  double baseRadius = params.radius;
  double topRadius = params.radius2;
  double height = params.length;
  discreteFunction[0][0] = 0.0;
  discreteFunction[0][1] = baseRadius;
  discreteFunction[1][0] = height;
  discreteFunction[1][1] = topRadius;

  axom::klee::Geometry vorGeometry =
    createGeometry_Vor(vorBase, vorDirection, discreteFunction);

  axom::klee::Shape vorShape("cone", "CONE", {}, {}, vorGeometry);

  return vorShape;
}

axom::klee::Shape createShape_Tet()
{
  axom::klee::TransformableGeometryProperties prop {
    axom::klee::Dimensions::Three,
    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  std::shared_ptr<axom::klee::Scale> scaleOp;
  if(!params.scaleFactors.empty())
  {
    scaleOp = std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                                  params.scaleFactors[1],
                                                  params.scaleFactors[2],
                                                  prop);
  }

  const double len = params.length;
  const double h = len / 3;  // To center the shape.
  const Point3D a {0.0 - h, 0.0 - h, 0.0 - h};
  const Point3D b {len - h, 0.0 - h, 0.0 - h};
  const Point3D c {0.0 - h, len - h, 0.0 - h};
  const Point3D d {0.0 - h, 0.0 - h, len - h};
  const primal::Tetrahedron<double, 3> tet {a, b, c, d};

  axom::klee::Geometry tetGeometry(prop, tet, scaleOp);
  axom::klee::Shape tetShape("tet", "TET", {}, {}, tetGeometry);

  return tetShape;
}

axom::klee::Shape createShape_Hex()
{
  axom::klee::TransformableGeometryProperties prop {
    axom::klee::Dimensions::Three,
    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  std::shared_ptr<axom::klee::Scale> scaleOp;
  if(!params.scaleFactors.empty())
  {
    scaleOp = std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                                  params.scaleFactors[1],
                                                  params.scaleFactors[2],
                                                  prop);
  }

  const double len = params.length;
  const double xhl = 0.5 * len * 0.8;
  const double yhl = 0.5 * len * 1.0;
  const double zhl = 0.5 * len * 1.2;
  // clang-format off
  const Point3D p {-xhl, -yhl, -zhl};
  const Point3D q { xhl, -yhl, -zhl};
  const Point3D r { xhl,  yhl, -zhl};
  const Point3D s {-xhl,  yhl, -zhl};
  const Point3D t {-xhl, -yhl,  zhl};
  const Point3D u { xhl, -yhl,  zhl};
  const Point3D v { xhl,  yhl,  zhl};
  const Point3D w {-xhl,  yhl,  zhl};
  // clang-format on
  const primal::Hexahedron<double, 3> hex {p, q, r, s, t, u, v, w};

  axom::klee::Geometry hexGeometry(prop, hex, scaleOp);
  axom::klee::Shape hexShape("hex", "HEX", {}, {}, hexGeometry);

  return hexShape;
}

axom::klee::Shape createShape_Plane()
{
  axom::klee::TransformableGeometryProperties prop {
    axom::klee::Dimensions::Three,
    axom::klee::LengthUnit::unspecified};

  SLIC_ASSERT(params.scaleFactors.empty() || params.scaleFactors.size() == 3);
  std::shared_ptr<axom::klee::Scale> scaleOp;
  if(!params.scaleFactors.empty())
  {
    scaleOp = std::make_shared<axom::klee::Scale>(params.scaleFactors[0],
                                                  params.scaleFactors[1],
                                                  params.scaleFactors[2],
                                                  prop);
  }

  // Create a plane crossing center of mesh.  No matter the normal,
  // it cuts the mesh in half.
  Point3D center {0.5 *
                  (primal::NumericArray<double, 3>(params.boxMins.data()) +
                   primal::NumericArray<double, 3>(params.boxMaxs.data()))};
  primal::Vector<double, 3> normal(params.direction.data());
  const primal::Plane<double, 3> plane {normal, center, true};

  axom::klee::Geometry planeGeometry(prop, plane, scaleOp);
  axom::klee::Shape planeShape("plane", "PLANE", {}, {}, planeGeometry);

  return planeShape;
}

//!@brief Create a ShapeSet with a single shape.
axom::klee::ShapeSet createShapeSet(const axom::klee::Shape& shape)
{
  axom::klee::ShapeSet shapeSet;
  shapeSet.setShapes(std::vector<axom::klee::Shape> {shape});
  shapeSet.setDimensions(axom::klee::Dimensions::Three);

  return shapeSet;
}

double volumeOfTetMesh(
  const axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& tetMesh)
{
  using TetType = axom::primal::Tetrahedron<double, 3>;
  if(0)
  {
    std::ofstream os("tets.js");
    tetMesh.getSidreGroup()->print(os);
  }
  axom::StackArray<axom::IndexType, 1> nodesShape {tetMesh.getNumberOfNodes()};
  axom::ArrayView<const double> x(tetMesh.getCoordinateArray(0), nodesShape);
  axom::ArrayView<const double> y(tetMesh.getCoordinateArray(1), nodesShape);
  axom::ArrayView<const double> z(tetMesh.getCoordinateArray(2), nodesShape);
  const axom::IndexType cellCount = tetMesh.getNumberOfCells();
  axom::Array<double> tetVolumes(cellCount, cellCount);
  double meshVolume = 0.0;
  for(axom::IndexType ic = 0; ic < cellCount; ++ic)
  {
    const axom::IndexType* nodeIds = tetMesh.getCellNodeIDs(ic);
    TetType tet;
    for(int j = 0; j < 4; ++j)
    {
      auto cornerNodeId = nodeIds[j];
      tet[j][0] = x[cornerNodeId];
      tet[j][1] = y[cornerNodeId];
      tet[j][2] = z[cornerNodeId];
    }
    meshVolume += tet.volume();
  }
  return meshVolume;
}

#if defined(AXOM_USE_MFEM)
/*!
  @brief Return the element volumes as a sidre::View.

  If it doesn't exist, allocate and compute it.
  \post The volume data is in \c dc->GetNamedBuffer(volFieldName).

  Most of this is lifted from IntersectionShaper::runShapeQueryImpl.
*/
template <typename ExecSpace>
axom::sidre::View* getElementVolumes(
  sidre::MFEMSidreDataCollection* dc,
  const std::string& volFieldName = std::string("elementVolumes"))
{
  using HexahedronType = axom::primal::Hexahedron<double, 3>;

  axom::sidre::View* volSidreView = dc->GetNamedBuffer(volFieldName);
  if(volSidreView == nullptr)
  {
    mfem::Mesh* mesh = dc->GetMesh();
    const axom::IndexType cellCount = mesh->GetNE();

    constexpr int NUM_VERTS_PER_HEX = 8;
    constexpr int NUM_COMPS_PER_VERT = 3;
    constexpr double ZERO_THRESHOLD = 1.e-10;

    axom::Array<Point3D> vertCoords(cellCount * NUM_VERTS_PER_HEX,
                                    cellCount * NUM_VERTS_PER_HEX);
    auto vertCoordsView = vertCoords.view();

    // This runs only only on host, because the mfem::Mesh only uses host memory, I think.
    for(axom::IndexType cellIdx = 0; cellIdx < cellCount; ++cellIdx)
    {
      // Get the indices of this element's vertices
      mfem::Array<int> verts;
      mesh->GetElementVertices(cellIdx, verts);
      SLIC_ASSERT(verts.Size() == NUM_VERTS_PER_HEX);

      // Get the coordinates for the vertices
      for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
      {
        int vertIdx = cellIdx * NUM_VERTS_PER_HEX + j;
        for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
        {
          vertCoordsView[vertIdx][k] = (mesh->GetVertex(verts[j]))[k];
        }
      }
    }

    // Set vertex coords to zero if within threshold.
    // (I don't know why we do this.  I'm following examples.)
    axom::ArrayView<double> flatCoordsView(
      (double*)vertCoords.data(),
      vertCoords.size() * Point3D::dimension());
    assert(flatCoordsView.size() == cellCount * NUM_VERTS_PER_HEX * 3);
    axom::for_all<ExecSpace>(
      cellCount * 3,
      AXOM_LAMBDA(axom::IndexType i) {
        if(axom::utilities::isNearlyEqual(flatCoordsView[i], 0.0, ZERO_THRESHOLD))
        {
          flatCoordsView[i] = 0.0;
        }
      });

    // Initialize hexahedral elements.
    axom::Array<HexahedronType> hexes(cellCount, cellCount);
    auto hexesView = hexes.view();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellIdx) {
        // Set each hexahedral element vertices
        hexesView[cellIdx] = HexahedronType();
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIndex = (cellIdx * NUM_VERTS_PER_HEX) + j;
          auto& hex = hexesView[cellIdx];
          hex[j] = vertCoordsView[vertIndex];
        }
      });  // end of loop to initialize hexahedral elements and bounding boxes

    // Allocate and populate cell volumes.
    volSidreView = dc->AllocNamedBuffer(volFieldName, cellCount);
    axom::ArrayView<double> volView(volSidreView->getData(),
                                    volSidreView->getNumElements());
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellIdx) {
        volView[cellIdx] = hexesView[cellIdx].volume();
      });
  }

  return volSidreView;
}
#endif

/*!
  @brief Return the element volumes as a sidre::View containing
  the volumes in an array.

  If it doesn't exist, allocate and compute it.
  \post The volume data is in the blueprint field \c volFieldName.

  Most of this is lifted from IntersectionShaper::runShapeQueryImpl.
*/
template <typename ExecSpace>
axom::sidre::View* getElementVolumes(
  sidre::Group* meshGrp,
  const std::string& volFieldName = std::string("elementVolumes"))
{
  using HexahedronType = axom::primal::Hexahedron<double, 3>;

  axom::sidre::View* volSidreView = nullptr;

  const auto fieldPath = axom::fmt::format("fields/{}", volFieldName);
  if(meshGrp->hasGroup(fieldPath))
  {
    sidre::Group* fieldGrp = meshGrp->getGroup(fieldPath);
    volSidreView = fieldGrp->getView("values");
  }
  else
  {
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(meshGrp,
                                                                topoName);

    const axom::IndexType cellCount = mesh.getNumberOfCells();

    constexpr int NUM_VERTS_PER_HEX = 8;
    constexpr int NUM_COMPS_PER_VERT = 3;
    constexpr double ZERO_THRESHOLD = 1.e-10;

    axom::Array<Point3D> vertCoords(cellCount * NUM_VERTS_PER_HEX,
                                    cellCount * NUM_VERTS_PER_HEX);
    auto vertCoordsView = vertCoords.view();

    // This runs only only on host, because the mfem::Mesh only uses host memory, I think.
    for(axom::IndexType cellIdx = 0; cellIdx < cellCount; ++cellIdx)
    {
      // Get the indices of this element's vertices
      axom::IndexType* verts = mesh.getCellNodeIDs(cellIdx);
      mesh.getCellNodeIDs(cellIdx, verts);

      // Get the coordinates for the vertices
      for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
      {
        int vertIdx = cellIdx * NUM_VERTS_PER_HEX + j;
        for(int k = 0; k < NUM_COMPS_PER_VERT; k++)
        {
          vertCoordsView[vertIdx][k] = mesh.getNodeCoordinate(verts[j], k);
        }
      }
    }

    // Set vertex coords to zero if within threshold.
    // (I don't know why we do this.  I'm following examples.)
    axom::ArrayView<double> flatCoordsView(
      (double*)vertCoords.data(),
      vertCoords.size() * Point3D::dimension());
    assert(flatCoordsView.size() == cellCount * NUM_VERTS_PER_HEX * 3);
    axom::for_all<ExecSpace>(
      cellCount * 3,
      AXOM_LAMBDA(axom::IndexType i) {
        if(axom::utilities::isNearlyEqual(flatCoordsView[i], 0.0, ZERO_THRESHOLD))
        {
          flatCoordsView[i] = 0.0;
        }
      });

    // Initialize hexahedral elements.
    axom::Array<HexahedronType> hexes(cellCount, cellCount);
    auto hexesView = hexes.view();
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellIdx) {
        // Set each hexahedral element vertices
        hexesView[cellIdx] = HexahedronType();
        for(int j = 0; j < NUM_VERTS_PER_HEX; ++j)
        {
          int vertIndex = (cellIdx * NUM_VERTS_PER_HEX) + j;
          auto& hex = hexesView[cellIdx];
          hex[j] = vertCoordsView[vertIndex];
        }
      });  // end of loop to initialize hexahedral elements and bounding boxes

    // Allocate and populate cell volumes.
    axom::sidre::Group* fieldGrp = meshGrp->createGroup(fieldPath);
    fieldGrp->createViewString("topology", topoName);
    fieldGrp->createViewString("association", "element");
    fieldGrp->createViewString("volume_dependent", "true");
    volSidreView =
      fieldGrp->createViewAndAllocate("values",
                                      axom::sidre::detail::SidreTT<double>::id,
                                      cellCount);
    axom::IndexType shape2d[] = {cellCount, 1};
    volSidreView->reshapeArray(2, shape2d);
    axom::ArrayView<double> volView(volSidreView->getData(),
                                    volSidreView->getNumElements());
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType cellIdx) {
        volView[cellIdx] = hexesView[cellIdx].volume();
      });
  }

  return volSidreView;
}

#if defined(AXOM_USE_MFEM)
/*!
  @brief Return global sum of volume of the given material.
*/
template <typename ExecSpace>
double sumMaterialVolumes(sidre::MFEMSidreDataCollection* dc,
                          const std::string& material)
{
  mfem::Mesh* mesh = dc->GetMesh();
  int const cellCount = mesh->GetNE();

  // Get cell volumes from dc.
  axom::sidre::View* elementVols = getElementVolumes<ExecSpace>(dc);
  axom::ArrayView<double> elementVolsView(elementVols->getData(),
                                          elementVols->getNumElements());

  // Get material volume fractions
  const std::string materialFieldName =
    axom::fmt::format("vol_frac_{}", material);
  mfem::GridFunction* volFracGf = dc->GetField(materialFieldName);
  axom::ArrayView<double> volFracGfArrayView(volFracGf->GetData(),
                                             volFracGf->Size());
  axom::quest::TempArrayView<ExecSpace> volFracView(volFracGfArrayView, true);

  using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<ReducePolicy, double> localVol(0);
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType i) {
      localVol += volFracView[i] * elementVolsView[i];
    });

  double globalVol = localVol.get();
  #ifdef AXOM_USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &globalVol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  return globalVol;
}
#endif

template <typename ExecSpace>
double sumMaterialVolumes(sidre::Group* meshGrp, const std::string& material)
{
  conduit::Node meshNode;
  meshGrp->createNativeLayout(meshNode);
#if defined(AXOM_DEBUG)
  conduit::Node info;
  conduit::blueprint::mesh::verify(meshNode, info);
  SLIC_ASSERT(conduit::blueprint::mesh::verify(meshNode, info));
#endif
  std::string topoPath = axom::fmt::format("topologies/{}", topoName);
  conduit::Node& topoNode = meshNode.fetch_existing(topoPath);
  const int cellCount = conduit::blueprint::mesh::topology::length(topoNode);

  // Get cell volumes from meshGrp.
  const std::string volsName = "vol_" + material;
  axom::sidre::View* elementVols =
    getElementVolumes<ExecSpace>(meshGrp, volsName);
  axom::ArrayView<double> elementVolsView(elementVols->getData(),
                                          elementVols->getNumElements());

  // Get material volume fractions
  const auto vfFieldName = axom::fmt::format("vol_frac_{}", material);
  const auto vfFieldValuesPath =
    axom::fmt::format("fields/{}/values", vfFieldName);
  axom::sidre::View* volFrac = meshGrp->getView(vfFieldValuesPath);
  axom::ArrayView<double> volFracArrayView(volFrac->getArray(), cellCount);
  axom::quest::TempArrayView<ExecSpace> volFracView(volFracArrayView, true);

  using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<ReducePolicy, double> localVol(0);
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType i) {
      localVol += volFracView[i] * elementVolsView[i];
    });

  double globalVol = localVol.get();
#ifdef AXOM_USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &globalVol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return globalVol;
}

/// Write blueprint mesh to disk
void saveMesh(const conduit::Node& mesh, const std::string& filename)
{
  AXOM_ANNOTATE_SCOPE("save mesh (conduit)");

#ifdef AXOM_USE_MPI
  conduit::relay::mpi::io::blueprint::save_mesh(mesh,
                                                filename,
                                                "hdf5",
                                                MPI_COMM_WORLD);
#else
  conduit::relay::io::blueprint::save_mesh(mesh, filename, "hdf5");
#endif
}

/// Write blueprint mesh to disk
void saveMesh(const sidre::Group& mesh, const std::string& filename)
{
  AXOM_ANNOTATE_SCOPE("save mesh (sidre)");

  conduit::Node tmpMesh;
  mesh.createNativeLayout(tmpMesh);
  {
    conduit::Node info;
#ifdef AXOM_USE_MPI
    if(!conduit::blueprint::mpi::verify("mesh", tmpMesh, info, MPI_COMM_WORLD))
#else
    if(!conduit::blueprint::verify("mesh", tmpMesh, info))
#endif
    {
      SLIC_INFO("Invalid blueprint for mesh: \n" << info.to_yaml());
      slic::flushStreams();
      assert(false);
    }
    // info.print();
  }
  saveMesh(tmpMesh, filename);
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::utilities::raii::MPIWrapper mpi_raii_wrapper(argc, argv);
  const int my_rank = mpi_raii_wrapper.my_rank();

  initializeLogger();

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  axom::CLI::App app {"Driver for Klee shaping query"};

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
    finalizeLogger();

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(retval);
  }

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
    params.annotationMode);

  AXOM_ANNOTATE_BEGIN("quest example for shaping primals");
  AXOM_ANNOTATE_BEGIN("init");

  // Storage for the shape geometry meshes.
  sidre::DataStore ds;

  //---------------------------------------------------------------------------
  // Create simple ShapeSet for the example.
  //---------------------------------------------------------------------------
  axom::klee::ShapeSet shapeSet;
  switch(params.boxDim)
  {
  case 2:
    shapeSet = create2DShapeSet(ds);
    break;
  case 3:
    if(params.testShape == "tetmesh")
    {
      shapeSet = createShapeSet(createShape_TetMesh(ds));
    }
    else if(params.testShape == "tet")
    {
      shapeSet = createShapeSet(createShape_Tet());
    }
    else if(params.testShape == "hex")
    {
      shapeSet = createShapeSet(createShape_Hex());
    }
    else if(params.testShape == "sphere")
    {
      shapeSet = createShapeSet(createShape_Sphere());
    }
    else if(params.testShape == "cyl")
    {
      shapeSet = createShapeSet(createShape_Cylinder());
    }
    else if(params.testShape == "cone")
    {
      shapeSet = createShapeSet(createShape_Cone());
    }
    else if(params.testShape == "vor")
    {
      shapeSet = createShapeSet(createShape_Vor());
    }
    else if(params.testShape == "plane")
    {
      shapeSet = createShapeSet(createShape_Plane());
    }
    break;
  }

  // Save the discrete shapes for viz and testing.
  auto* shapeMeshGroup = ds.getRoot()->createGroup("shapeMeshGroup");
  std::vector<std::shared_ptr<axom::mint::Mesh>> discreteShapeMeshes;
  for(const auto& shape : shapeSet.getShapes())
  {
    axom::quest::DiscreteShape dShape(shape, shapeMeshGroup);
    auto dMesh =
      std::dynamic_pointer_cast<axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>>(
        dShape.createMeshRepresentation());
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format("Shape '{}' discrete geometry has {} cells",
                        shape.getName(),
                        dMesh->getNumberOfCells())));

    discreteShapeMeshes.push_back(dMesh);

    if(!params.outputFile.empty())
    {
      std::string shapeFileName = params.outputFile + ".shape";
      conduit::Node tmpNode, info;
      dMesh->getSidreGroup()->createNativeLayout(tmpNode);
      conduit::relay::io::blueprint::save_mesh(tmpNode, shapeFileName, "hdf5");
    }
  }

  const klee::Dimensions shapeDim = shapeSet.getDimensions();

  // Apply error checking
#ifndef AXOM_USE_C2C
  SLIC_ERROR_IF(shapeDim == klee::Dimensions::Two,
                "Shaping with contour files requires an Axom configured with "
                "the C2C library");
#endif

  axom::IndexType cellCount = -1;

#if defined(AXOM_USE_MFEM)
  std::shared_ptr<sidre::MFEMSidreDataCollection> shapingDC;
  if(params.useMfem())
  {
    AXOM_ANNOTATE_BEGIN("load mesh");
    //---------------------------------------------------------------------------
    // Load the computational mesh
    // originalMeshDC is the input MFEM mesh
    // It's converted to shapingDC below, and that is used for shaping calls.
    //---------------------------------------------------------------------------
    auto originalMeshDC = params.loadComputationalMesh();

    //---------------------------------------------------------------------------
    // Set up DataCollection for shaping
    // shapingDC is a "copy" of originalMeshDC, the MFEM mesh.
    // It's created empty then populated by the SetMesh call.
    // shapingMesh and parallelMesh are some kind of temporary versions of originalMeshDC.
    //---------------------------------------------------------------------------
    mfem::Mesh* shapingMesh = nullptr;
    constexpr bool dc_owns_data = true;
    shapingDC = std::make_shared<sidre::MFEMSidreDataCollection>("shaping",
                                                                 shapingMesh,
                                                                 dc_owns_data);
    {
      shapingDC->SetMeshNodesName("positions");

      // With MPI, loadComputationalMesh returns a parallel mesh.
      mfem::ParMesh* parallelMesh =
        dynamic_cast<mfem::ParMesh*>(originalMeshDC->GetMesh());
      shapingMesh = (parallelMesh != nullptr)
        ? new mfem::ParMesh(*parallelMesh)
        : new mfem::Mesh(*originalMeshDC->GetMesh());
      shapingDC->SetMesh(shapingMesh);
    }
    AXOM_ANNOTATE_END("load mesh");
    printMfemMeshInfo(shapingDC->GetMesh(), "Loaded MFEM mesh");

    cellCount = shapingMesh->GetNE();
  }
#else
  SLIC_ERROR_IF(params.useMfem(),
                "Cannot use MFEM mesh due to Axom configuration.  Please use "
                "Blueprint mesh or configure with MFEM and "
                "-DAXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION.");
#endif

  axom::sidre::Group* compMeshGrp = nullptr;
  if(params.useBlueprint())
  {
    compMeshGrp = createBoxMesh(ds.getRoot()->createGroup("compMesh"));
    conduit::Node meshNode;
    compMeshGrp->createNativeLayout(meshNode);
    SLIC_INFO(axom::fmt::format("{:-^80}", "Generated Blueprint mesh"));
    meshNode.print();
    const conduit::Node& topoNode = meshNode["topologies"][topoName];
    cellCount = conduit::blueprint::mesh::topology::length(topoNode);
  }

  // TODO Port to GPUs.  Shaper should be data-parallel, but data may not be on devices yet.
  using ExecSpace = typename axom::SEQ_EXEC;
  using ReducePolicy = typename axom::execution_space<ExecSpace>::reduce_policy;

  //---------------------------------------------------------------------------
  // Initialize the shaping query object
  //---------------------------------------------------------------------------
  AXOM_ANNOTATE_BEGIN("setup shaping problem");
  std::shared_ptr<quest::IntersectionShaper> shaper = nullptr;
  if(params.useBlueprint())
  {
    shaper = std::make_shared<quest::IntersectionShaper>(shapeSet, compMeshGrp);
  }
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    shaper =
      std::make_shared<quest::IntersectionShaper>(shapeSet, shapingDC.get());
  }
#endif
  SLIC_ASSERT(shaper != nullptr);

  // Set generic parameters for the base Shaper instance
  shaper->setVertexWeldThreshold(params.weldThresh);
  shaper->setVerbosity(params.isVerbose());
  if(params.percentError > 0.)
  {
    shaper->setPercentError(params.percentError);
    shaper->setRefinementType(quest::DiscreteShape::RefinementDynamic);
  }

  // Associate any fields that begin with "vol_frac" with "material" so when
  // the data collection is written, a matset will be created.
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    shaper->getDC()->AssociateMaterialSet("vol_frac", "material");
  }
#endif

  // Set specific parameters here for IntersectionShaper
  shaper->setLevel(params.refinementLevel);
  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Setting IntersectionShaper policy to '{}'",
                      axom::runtime_policy::policyToName(params.policy))));
  shaper->setExecPolicy(params.policy);

  if(!params.backgroundMaterial.empty())
  {
    shaper->setFreeMaterialName(params.backgroundMaterial);
  }

  AXOM_ANNOTATE_END("setup shaping problem");
  AXOM_ANNOTATE_END("init");

  //---------------------------------------------------------------------------
  // Process each of the shapes
  //---------------------------------------------------------------------------
  SLIC_INFO(axom::fmt::format("{:=^80}", "Shaping loop"));
  AXOM_ANNOTATE_BEGIN("shaping");
  for(const auto& shape : shapeSet.getShapes())
  {
    const std::string shapeFormat = shape.getGeometry().getFormat();
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format("Processing shape '{}' of material '{}' (format '{}')",
                        shape.getName(),
                        shape.getMaterial(),
                        shapeFormat)));

    // Load the shape from file. This also applies any transformations.
    shaper->loadShape(shape);
    slic::flushStreams();

    // Generate a spatial index over the shape
    shaper->prepareShapeQuery(shapeDim, shape);
    slic::flushStreams();

    // Query the mesh against this shape
    shaper->runShapeQuery(shape);
    slic::flushStreams();

    // Apply the replacement rules for this shape against the existing materials
    shaper->applyReplacementRules(shape);
    slic::flushStreams();

    // Finalize data structures associated with this shape and spatial index
    shaper->finalizeShapeQuery();
    slic::flushStreams();
  }
  AXOM_ANNOTATE_END("shaping");

  //---------------------------------------------------------------------------
  // After shaping in all shapes, generate/adjust the material volume fractions
  //---------------------------------------------------------------------------
  AXOM_ANNOTATE_BEGIN("adjust");
  SLIC_INFO(
    axom::fmt::format("{:=^80}",
                      "Generating volume fraction fields for materials"));

  shaper->adjustVolumeFractions();

  //---------------------------------------------------------------------------
  // Compute and print volumes of each material's volume fraction
  //---------------------------------------------------------------------------
  using axom::utilities::string::startsWith;
  if(params.useBlueprint())
  {
    std::vector<std::string> materialNames = shaper->getMaterialNames();
    for(const auto& materialName : materialNames)
    {
      // Compute and print volume of material.
      const double volume =
        sumMaterialVolumes<axom::SEQ_EXEC>(compMeshGrp, materialName);
      SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                  "Volume of material '{}' is {:.6Lf}",
                                  materialName,
                                  volume));
    }
  }
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    for(auto& kv : shaper->getDC()->GetFieldMap())
    {
      if(startsWith(kv.first, "vol_frac_"))
      {
        const auto mat_name = kv.first.substr(9);
        auto* gf = kv.second;

        mfem::ConstantCoefficient one(1.0);
        mfem::LinearForm vol_form(gf->FESpace());
        vol_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
        vol_form.Assemble();

        const double volume = shaper->allReduceSum(*gf * vol_form);

        SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                    "Volume of material '{}' is {:.6Lf}",
                                    mat_name,
                                    volume));
      }
    }
  }
#endif
  AXOM_ANNOTATE_END("adjust");

  int failCounts = 0;

  axom::sidre::Group* volFracGroups = nullptr;
  if(params.useBlueprint())
  {
    volFracGroups = compMeshGrp->getGroup("matsets/material/volume_fractions");
  }
#if defined(AXOM_USE_MFEM)
  if(params.useMfem())
  {
    volFracGroups =
      shapingDC->GetBPGroup()->getGroup("matsets/material/volume_fractions");
  }
#endif

  //---------------------------------------------------------------------------
  // Correctness test: volume fractions should be in [0,1].
  //---------------------------------------------------------------------------
  RAJA::ReduceSum<ReducePolicy, axom::IndexType> rangeViolationCount(0);
  for(axom::sidre::Group& materialGroup : volFracGroups->groups())
  {
    axom::sidre::View* values = materialGroup.getView("value");
    double* volFracData = values->getArray();
    axom::ArrayView<double> volFracDataView(volFracData, cellCount);
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) {
        bool bad = volFracDataView[i] < 0.0 || volFracDataView[i] > 1.0;
        rangeViolationCount += bad;
      });
  }

  failCounts += (rangeViolationCount.get() != 0);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Count of volume fractions outside of [0,1]: {}.",
                      rangeViolationCount.get())));
  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Correctness test: volume fractions in each cell should sum to 1.0.
  //---------------------------------------------------------------------------
  axom::Array<double> volSums(cellCount);
  volSums.fill(0.0);
  axom::ArrayView<double> volSumsView = volSums.view();
  for(axom::sidre::Group& materialGroup : volFracGroups->groups())
  {
    axom::sidre::View* values = materialGroup.getView("value");
    double* volFracData = values->getArray();
    axom::ArrayView<double> volFracDataView(volFracData, cellCount);
    axom::for_all<ExecSpace>(
      cellCount,
      AXOM_LAMBDA(axom::IndexType i) { volSumsView[i] += volFracDataView[i]; });
  }
  RAJA::ReduceSum<ReducePolicy, axom::IndexType> nonUnitSums(0);
  axom::for_all<ExecSpace>(
    cellCount,
    AXOM_LAMBDA(axom::IndexType i) {
      bool bad = !axom::utilities::isNearlyEqual(volSums[i], 0.0);
      nonUnitSums += bad;
    });

  failCounts += (nonUnitSums.get() != 0);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Count non-unit volume fraction sums: {}.",
                      nonUnitSums.get())));
  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Correctness test: shape volume in shapingMesh should match volume of the
  // shape mesh for closes shape.
  //---------------------------------------------------------------------------
  auto* meshVerificationGroup = ds.getRoot()->createGroup("meshVerification");
  for(const auto& shape : shapeSet.getShapes())
  {
    axom::quest::DiscreteShape dShape(shape, meshVerificationGroup);
    auto shapeMesh =
      std::dynamic_pointer_cast<axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>>(
        dShape.createMeshRepresentation());
    double shapeMeshVol = volumeOfTetMesh(*shapeMesh);
    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format("Shape '{}' discrete geometry has {} cells",
                        shape.getName(),
                        shapeMesh->getNumberOfCells())));

    const std::string& materialName = shape.getMaterial();
    double shapeVol = -1;
    if(params.useBlueprint())
    {
      shapeVol = sumMaterialVolumes<axom::SEQ_EXEC>(compMeshGrp, materialName);
    }
#if defined(AXOM_USE_MFEM)
    if(params.useMfem())
    {
      shapeVol =
        sumMaterialVolumes<axom::SEQ_EXEC>(shapingDC.get(), materialName);
    }
#endif
    double correctShapeVol =
      params.testShape == "plane" ? params.boxMeshVolume() / 2 : shapeMeshVol;
    SLIC_ASSERT(correctShapeVol > 0.0);  // Indicates error in the test setup.
    double diff = shapeVol - correctShapeVol;

    bool err = !axom::utilities::isNearlyEqualRelative(shapeVol,
                                                       correctShapeVol,
                                                       1e-6,
                                                       1e-8);
    failCounts += err;

    SLIC_INFO(axom::fmt::format(
      "{:-^80}",
      axom::fmt::format(
        "Material '{}' in shape '{}' has volume {} vs {}, diff of {}, {}.",
        materialName,
        shape.getName(),
        shapeVol,
        correctShapeVol,
        diff,
        (err ? "ERROR" : "OK"))));
  }
  slic::flushStreams();
  ds.getRoot()->destroyGroupAndData("meshVerification");

  //---------------------------------------------------------------------------
  // Save meshes and fields
  //---------------------------------------------------------------------------

  if(!params.outputFile.empty())
  {
    std::string fileName = params.outputFile + ".volfracs";
    if(params.useBlueprint())
    {
      saveMesh(*compMeshGrp, fileName);
      SLIC_INFO(axom::fmt::format("{:-^80}", "Wrote output mesh " + fileName));
    }
#if defined(AXOM_USE_MFEM)
    else
    {
      shaper->getDC()->Save(fileName, sidre::Group::getDefaultIOProtocol());
    }
#endif
  }

  shaper.reset();

  //---------------------------------------------------------------------------
  // Cleanup and exit
  //---------------------------------------------------------------------------
  SLIC_INFO(axom::fmt::format("{:-^80}", ""));
  slic::flushStreams();

  AXOM_ANNOTATE_END("quest shaping example");

  finalizeLogger();

  return failCounts;
}
