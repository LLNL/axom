// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_signed_distance_mfem.cpp
 * \brief Computes signed distance field over an mfem mesh
 */

/*
 *  TODO:
 *  [] Add mfem mesh to data collection
 *  [] Add distance field as H1 grid function
 *  [] Call SignedDistance on all DOFs of distance field
 *  [] Output results
 *  [] Add ability to transform mesh (scale and translate) 
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/spin.hpp"
#include "axom/mint.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include "fmt/fmt.hpp"
#include "CLI11/CLI11.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#else
using MPI_Comm = int;
#endif

#ifndef AXOM_USE_MFEM
  #error "This example depends on mfem"
#endif
#include "mfem.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>
#include <iomanip>  // for setprecision()

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace slic = axom::slic;

/** Struct to parse and store the input parameters */
struct Input
{
public:
  std::string stlFile;
  std::string mfemFile;
  int distanceOrder {2};

  mfem::Mesh* mesh {nullptr};
  mfem::DataCollection* dc {nullptr};

private:
  bool m_verboseOutput {false};
  bool m_use_batched_query {false};

public:
  Input() = default;
  ~Input()
  {
    delete dc;
    dc = nullptr;

    delete mesh;
    mesh = nullptr;

    quest::signed_distance_finalize();
  }

  bool isVerbose() const { return m_verboseOutput; }
  bool useBatchedQuery() const { return m_use_batched_query; }

  void loadMfemMesh()
  {
    SLIC_INFO(fmt::format("Loading mfem mesh '{}'", mfemFile));

    const int generateEdges = 1;
    const int refinement = 1;
    mesh = new mfem::Mesh(mfemFile.c_str(), generateEdges, refinement);
    mesh->EnsureNodes();
  }

  void setupSignedDistance()
  {
    SLIC_INFO(fmt::format("Initializing signed distance structures over '{}'",
                          stlFile));

#ifdef AXOM_USE_MPI
    MPI_Comm comm {MPI_COMM_WORLD};
#else
    MPI_Comm comm {0};
#endif

    //quest::signed_distance_set_max_levels(level);
    //quest::signed_distance_set_max_occupancy(occupancy);
    //quest::signed_distance_set_compute_signs(shouldCompute);
    quest::signed_distance_init(stlFile, comm);
  }

  template <int DIM>
  void printMFEMMeshBoundingBox()
  {
    using BBox = primal::BoundingBox<double, DIM>;
    using Pt = typename BBox::PointType;

    mfem::Vector bmin, bmax;
    mesh->GetBoundingBox(bmin, bmax);
    SLIC_INFO("Bounding box for mfem mesh: " << BBox(Pt(bmin.GetData(), DIM),
                                                     Pt(bmax.GetData(), DIM)));
  }

  template <int DIM>
  void printStlMeshBoundingBox()
  {
    using BBox = primal::BoundingBox<double, DIM>;
    using Pt = typename BBox::PointType;

    Pt bmin, bmax;
    quest::signed_distance_get_mesh_bounds(bmin.data(), bmax.data());
    SLIC_INFO("Bounding box for STL mesh: " << BBox(bmin, bmax));
  }

  void printStats()
  {
    const int dim = mesh->Dimension();

    switch(dim)
    {
    case 2:
      printMFEMMeshBoundingBox<2>();
      printStlMeshBoundingBox<2>();
      break;
    case 3:
      printMFEMMeshBoundingBox<3>();
      printStlMeshBoundingBox<3>();
      break;
    }

    {
      std::stringstream sstr;
      mesh->PrintInfo(sstr);
      SLIC_INFO("MFEM mesh info: \n" << sstr.str());
    }
  }

  void parse(int argc, char** argv, CLI::App& app)
  {
    app.add_option("-m, --mfemMesh", mfemFile, "Path to input mfem mesh")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-s, --stlFile", stlFile, "Path to surface mesh")
      ->required()
      ->check(CLI::ExistingFile);
    // maybe add an output filename?

    app
      .add_flag("-v,--verbose", m_verboseOutput, "Enable/disable verbose output")
      ->capture_default_str();

    app
      .add_flag("-k,--order",
                distanceOrder,
                "Polynomial order for computed distance field")
      ->capture_default_str()
      ->check(CLI::PositiveNumber);

    app.get_formatter()->column_width(35);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);
  }
};

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  int mpirank {0};
  int numranks {1};

#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
#endif

  axom::slic::UnitTestLogger logger;  // create & initialize logger
  // slic::debug::checksAreErrors = true;

  // Set up and parse command line arguments
  Input params;
  CLI::App app {
    "Computes the signed distance field to a surface over an mfem mesh"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const CLI::ParseError& e)
  {
    int retval = -1;
    if(mpirank == 0)
    {
      retval = app.exit(e);
    }

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(retval);
  }

  // load mfem mesh and signed distance data structures
  params.loadMfemMesh();
  params.setupSignedDistance();
  params.printStats();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
