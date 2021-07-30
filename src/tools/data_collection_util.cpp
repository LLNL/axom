// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"
#include "axom/primal.hpp"

#include "fmt/fmt.hpp"
#include "CLI11/CLI11.hpp"

// MFEM is required
#ifndef AXOM_USE_MFEM
  #error This example requires an Axom configured with MFEM
#endif
#include "mfem.hpp"

// MPI is optional
#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

namespace slic = axom::slic;
namespace sidre = axom::sidre;
namespace primal = axom::primal;

//------------------------------------------------------------------------------

/// Struct to help choose if we're loading a box (Cartesian) mesh or a file
enum class MeshForm : int
{
  Box,
  File
};

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string dcName {"mesh"};

  int uniformRefinements {0};
  int polynomialOrder {2};

  MeshForm meshForm {MeshForm::Box};

  // Box mesh parameters
  std::vector<double> boxMins;
  std::vector<double> boxMaxs;
  std::vector<int> boxResolution;
  int boxDim {-1};

  std::string mfemFile;

private:
  bool m_verboseOutput {false};

  /**
   * Fixes up parameters associated with a box mesh based on user-provided input
   */
  void fixBoxParams()
  {
    const bool dimProvided = boxDim > 0;
    const bool rangeProvided = !(boxMins.empty() && boxMaxs.empty());
    const bool resProvided = !boxResolution.empty();

    // First check: we need at least one of: dimension, mins/maxs, resolution
    if(!(dimProvided || rangeProvided || resProvided))
    {
      throw "Box mesh must have at least one of: dimension, mins and maxs or "
        "resolution to determine the dimension";
    }

    int szMins = boxMins.size();
    int szMaxs = boxMaxs.size();
    int szRes = boxResolution.size();

    // Error checking on provided dimension
    if(dimProvided)
    {
      if(boxDim < 1 || boxDim > 3)
      {
        throw "Dimension cannot exceed 3";
      }

      // Ensure that range and resolution have the right number of entries
      if(rangeProvided && (szMins != boxDim || szMaxs != boxDim))
      {
        throw "Bounding box has different dimension than provided dimension";
      }
      if(resProvided && (szRes != boxDim))
      {
        throw "Box resolution has different dimension than provided dimension";
      }
    }

    // Error checking on provided range; note: accounts for previous checks on dim
    if(rangeProvided && (szMins != szMaxs))
    {
      throw "Bounding box mins and maxs has different dimensions";
    }

    // Error checking on provided range and resolution; note: accounts for previous checks on dim and range
    if(rangeProvided && resProvided && (szMins != szRes))
    {
      throw "Bounding box mins and maxs has different dimensions than resolution";
    }

    // Now that we've checked errors, let's fill in missing data
    const int dim = dimProvided ? boxDim : (resProvided ? szRes : szMins);
    if(!dimProvided)
    {
      boxDim = dim;
    }

    if(!rangeProvided)  // set range to unit box
    {
      boxMins.resize(dim);
      boxMaxs.resize(dim);
      for(int i = 0; i < dim; ++i)
      {
        boxMins[i] = 0.;
        boxMaxs[i] = 1.;
      }
    }

    if(!resProvided)  // set resolution to 1
    {
      boxResolution.resize(dim);
      for(int i = 0; i < dim; ++i)
      {
        boxResolution[i] = 1;
      }
    }
  }

public:
  Input() = default;

  bool isVerbose() const { return m_verboseOutput; }

  void parse(int argc, char** argv, CLI::App& app)
  {
    app.add_option("-m, --mfem-file", mfemFile)
      ->description("Path to a mesh file in the mfem format")
      ->check(CLI::ExistingFile);

    app.add_flag("-v,--verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("-p,--polynomial-order", polynomialOrder)
      ->description("polynomial order of the generated mesh")
      ->capture_default_str()
      ->check(CLI::NonNegativeNumber);

    std::map<std::string, MeshForm> meshFormMap {{"box", MeshForm::Box},
                                                 {"file", MeshForm::File}};
    app.add_option("-f,--mesh-form", meshForm)
      ->description("Determines the input type -- either box or file")
      ->capture_default_str()
      ->transform(CLI::CheckedTransformer(meshFormMap, CLI::ignore_case));

    app.add_option("-l,--ref-level", uniformRefinements)
      ->description(
        "The number of uniform refinement levels to apply to the mesh")
      ->capture_default_str()
      ->check(CLI::NonNegativeNumber);

    // Optional bounding box for query region
    auto* minbb = app.add_option("--min", boxMins)
                    ->description("Min bounds for box mesh (x,y[,z])")
                    ->expected(2, 3);
    auto* maxbb = app.add_option("--max", boxMaxs)
                    ->description("Max bounds for box mesh (x,y[,z])")
                    ->expected(2, 3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    app.add_option("--res", boxResolution)
      ->description("Reesolution of the box mesh (i,j[,k])")
      ->expected(2, 3);

    app.add_option("-d,--dimension", boxDim)
      ->description("dimension of the box mesh")
      ->capture_default_str()
      ->check(CLI::PositiveNumber);

    app.get_formatter()->column_width(50);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);

    if(meshForm == MeshForm::Box)
    {
      fixBoxParams();  // Note: throws on error conditions
    }
    else
    {
      if(mfemFile.empty())
      {
        throw "'mfemFile' required when `--mesh-form == File`";
      }
    }
  }
};

mfem::Mesh* createBoxMesh(Input& params)
{
  // Create a background mesh
  // Generate an mfem Cartesian mesh, scaled to the bounding box range

  const int dim = params.boxDim;
  auto& lo = params.boxMins;
  auto& hi = params.boxMaxs;
  auto& res = params.boxResolution;

  mfem::Mesh* mesh = nullptr;

  // TODO: Convert to MakeCartesian2D and MakeCartesian3D when upgrading to mfem@4.3
  switch(dim)
  {
  case 2:
    SLIC_INFO(fmt::format(
      "Creating a box mesh of resolution {} and bounding box {}",
      primal::Point<int, 2>(res.data()),
      primal::BoundingBox<double, 2>(primal::Point<double, 2>(lo.data()),
                                     primal::Point<double, 2>(hi.data()))));

    mesh = new mfem::Mesh(res[0],
                          res[1],
                          mfem::Element::QUADRILATERAL,
                          false,
                          hi[0] - lo[0],
                          hi[1] - lo[1]);
    break;
  case 3:
    SLIC_INFO(fmt::format(
      "Creating a box mesh of resolution {} and bounding box {}",
      primal::Point<int, 3>(res.data()),
      primal::BoundingBox<double, 3>(primal::Point<double, 3>(lo.data()),
                                     primal::Point<double, 3>(hi.data()))));

    mesh = new mfem::Mesh(res[0],
                          res[1],
                          res[2],
                          mfem::Element::HEXAHEDRON,
                          false,
                          hi[0] - lo[0],
                          hi[1] - lo[1],
                          hi[2] - lo[2]);
    break;
  default:
    SLIC_ERROR("Only 2D and 3D meshes are currently supported.");
    break;
  }

  // Offset to the mesh to lie w/in the bounding box
  for(int i = 0; i < mesh->GetNV(); ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < dim; ++d)
    {
      v[d] += lo[d];
    }
  }

  // Ensure that mesh has high order nodes
  mesh->SetCurvature(params.polynomialOrder);

  return mesh;
}

/*!
 * \brief Utility function to initialize the logger
 */
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
  const int RLIMIT = 8;
  logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
  logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
#else
  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif  // AXOM_USE_MPI

  slic::addStreamToAllMsgLevels(logStream);
}

/*!
 * \brief Utility function to finalize the logger
 */
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
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

  initializeLogger();

  // Set up and parse command line arguments
  Input params;
  CLI::App app {"Utility tool to create a blueprint compliant data store"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const CLI::ParseError& e)
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

  // Create the mfem sidre data collection
  sidre::MFEMSidreDataCollection dc(params.dcName, nullptr, true);

  // Create or load the serial mfem mesh
  mfem::Mesh* mesh = nullptr;
  switch(params.meshForm)
  {
  case MeshForm::Box:
    mesh = createBoxMesh(params);
    break;

  case MeshForm::File:
    mesh = new mfem::Mesh(params.mfemFile.c_str(), 1, 1);
    break;
  }

  // Apply uniform refinements
  for(int i = 0; i < params.uniformRefinements; ++i)
  {
    mesh->UniformRefinement();
  }

  // Handle conversion to parallel mfem mesh
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  {
    int* partitioning = nullptr;
    int part_method = 0;
    mfem::Mesh* pmesh =
      new mfem::ParMesh(MPI_COMM_WORLD, *mesh, partitioning, part_method);
    delete[] partitioning;
    delete mesh;
    mesh = pmesh;
  }
#endif

  // Add mfem mesh to data collection
  dc.SetMeshNodesName("positions");
#ifdef MFEM_USE_MPI
  dc.SetMesh(MPI_COMM_WORLD, mesh);
#else
  dc.SetMesh(mesh);
#endif

  // TODO: Optionally convert to low order mesh ?

  SLIC_INFO(fmt::format("Saving mesh file '{}'", dc.GetCollectionName()));
#ifdef MFEM_USE_MPI
  dc.Save();
#endif

  // Cleanup and exit
  finalizeLogger();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
