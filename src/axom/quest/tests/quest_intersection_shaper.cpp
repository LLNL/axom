// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file
 * \brief Unit tests for quest's IntersectionShaper class replacement rules.
 *
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/klee.hpp"
#include "axom/quest.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"
// _quest_intersection_shaper_include_start
#include "axom/quest/IntersectionShaper.hpp"

#ifndef AXOM_USE_MFEM
  #error "Quest's IntersectionShaper tests on mfem meshes require mfem library."
#endif
#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif
#include <cmath>
#include <string>
#include <vector>
// _quest_intersection_shaper_include_end

//#define AXOM_SEQ_INTERSECTION_SHAPER_WORKS
#define GENERATE_BASELINES
#define VISUALIZE_DATASETS

std::vector<std::string> case1{
  "shaping/case1/case1_012.yaml",
  "shaping/case1/case1_021.yaml",
  "shaping/case1/case1_102.yaml",
  "shaping/case1/case1_120.yaml",
  "shaping/case1/case1_201.yaml",
  "shaping/case1/case1_210.yaml"
};

std::vector<std::string> case2{
  "shaping/case2/case2_012.yaml",
  "shaping/case2/case2_021.yaml",
  "shaping/case2/case2_102.yaml",
  "shaping/case2/case2_120.yaml",
  "shaping/case2/case2_201.yaml",
  "shaping/case2/case2_210.yaml"
};

std::vector<std::string> case3{
  "shaping/case3/case3_012.yaml",
  "shaping/case3/case3_021.yaml",
  "shaping/case3/case3_102.yaml",
  "shaping/case3/case3_120.yaml",
  "shaping/case3/case3_201.yaml",
  "shaping/case3/case3_210.yaml"
};

std::vector<std::string> case4{
  "shaping/case4/case4.yaml",
  "shaping/case4/case4_overwrite.yaml"
};

std::vector<std::string> case5{
  "shaping/case5.yaml"
};

using std::cout;
using std::endl;

namespace quest = axom::quest;
namespace slic = axom::slic;
namespace sidre = axom::sidre;
namespace klee = axom::klee;

std::string
pjoin(const std::string &path, const std::string &filename)
{
  return path + "/" + filename;
}

void
psplit(const std::string &filepath, std::string &path, std::string &filename)
{
  auto idx = filepath.rfind("/");
  if(idx != std::string::npos)
  {
    path = filepath.substr(0, idx);
    filename = filepath.substr(idx+1, filepath.size()-1-idx-1);
  }
  else
  {
    path = "";
    filename = filepath;
  }
}

std::string
dataDirectory()
{
// TEMP
  return "/usr/WS2/whitlocb/Axom/axom_data";

  return AXOM_DATA_DIR;
}

std::string
testData(const std::string &filename)
{
  return pjoin(dataDirectory(), filename);
}

std::string
baselineDir()
{
  return pjoin(pjoin(pjoin(dataDirectory(), "quest"), "regression"), "quest_intersection_shaper");
}

std::string
yaml_root(const std::string &filepath)
{
  std::string retval, path, filename;
  psplit(filepath, path, filename);
  auto idx = filename.rfind(".");
  if(idx != std::string::npos)
    retval = filename.substr(0, idx);
  else
    retval = filename;
  return retval;
}

mfem::GridFunction *
newGridFunction(mfem::Mesh *mesh)
{
  const int vfOrder = 0;
  const int dim = mesh->Dimension();
  mfem::L2_FECollection* coll =
    new mfem::L2_FECollection(vfOrder, dim, mfem::BasisType::Positive);
  mfem::FiniteElementSpace* fes = new mfem::FiniteElementSpace(mesh, coll);
  mfem::GridFunction* gf = new mfem::GridFunction(fes);
  gf->MakeOwner(coll);
  memset(gf->GetData(), 0, gf->Size() * sizeof(double));
  return gf;
}

void
makeTestMesh(sidre::MFEMSidreDataCollection &dc)
{
  int polynomialOrder = 1;
  double lo[] = {0., 0., -0.25};
  double hi[] = {1., 1., 0.};
  int celldims[] = {20, 20, 1};
  auto mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian3D(celldims[0],
                                                         celldims[1],
                                                         celldims[2],
                                                         mfem::Element::HEXAHEDRON,
                                                         hi[0] - lo[0],
                                                         hi[1] - lo[1],
                                                         hi[2] - lo[2],
                                                         false));
  mesh->SetCurvature(polynomialOrder);
  dc.SetMeshNodesName("positions");
  dc.SetMesh(mesh);

#if 0
  // This mode will make 2 clean materials in the mesh.
  mfem::GridFunction *mata = newGridFunction(mesh);
  mfem::GridFunction *matb = newGridFunction(mesh);
  for(int k = 0; k < celldims[2]; k++)
  {
    for(int j = 0; j < celldims[1]; j++)
    {
      for(int i = 0; i < celldims[0]; j++)
      {
        int id = k*celldims[1]*celldims[0] + j*celldims[0] + i;
        (*mata)(id) = (j < celldims[1]/2);
        (*matb)(id) = (j >= celldims[1]/2);
      }
    }
  }
  dc.RegisterField("vol_frac_a", mata);
  dc.RegisterField("vol_frac_b", matb);
#endif
}

bool
compareGF(const mfem::GridFunction *baselineGF, const mfem::GridFunction *newGF,
  double tolerance)
{
  bool retval = false;
  if(baselineGF->Size() == newGF->Size())
  {
    const double *bgf = baselineGF->GetData();
    const double *ngf = newGF->GetData();
    double maxdiff = 0.;
    for(int i = 0; i < baselineGF->Size(); i++)
    {
      double diff = fabs(bgf[i] - ngf[i]);
      maxdiff = std::max(maxdiff, diff);
      if(diff > tolerance)
      {
        cout << "Grid functions differ at index " << i << ". maxdiff=" << maxdiff << endl;
        retval = false;
        break;
      }
    }
    retval = true;
  }
  return retval;
}

// compare the volume fraction fields in the DCs.
bool
compareDC(sidre::MFEMSidreDataCollection &baselineDC,
  sidre::MFEMSidreDataCollection &newDC, double tolerance)
{
  for(auto it : baselineDC.GetFieldMap())
  {
    // Just compare vol_frac_ grid functions.
    if(it.first.find("vol_frac_") != std::string::npos)
    {
      // The field is a shape or volfrac.
      auto newGF = newDC.GetField(it.first);
      if(newGF)
      {
        if(!compareGF(it.second, newGF, tolerance))
        {
          cout << "Grid function " << it.first << " does not match baseline." << endl;
          return false;
        }
      }
      else
      {
        cout << "Grid function " << it.first << " does not exist in new data collection." << endl;
        return false;
      }
    }
  }
  return true;
}

// Save Sidre as VisIt
void
saveVisIt(const std::string &path, const std::string &filename, sidre::MFEMSidreDataCollection &dc)
{
  // Wrap mesh and grid functions in a VisItDataCollection and save it.
  mfem::VisItDataCollection vdc(filename, dc.GetMesh());
  if(!path.empty())
    vdc.SetPrefixPath(path);
  vdc.SetOwnData(false);
  vdc.SetFormat(mfem::DataCollection::SERIAL_FORMAT);
  for(auto it : dc.GetFieldMap())
  {
    if(it.first.find("vol_frac_") != std::string::npos)
      vdc.RegisterField(it.first, it.second);
  } 
  vdc.Save();
}

// Load VisIt as Sidre
void
loadVisIt(mfem::VisItDataCollection &vdc, sidre::MFEMSidreDataCollection &dc)
{
  // Wrap mesh and grid functions in a VisItDataCollection and save it.
  vdc.SetFormat(mfem::DataCollection::SERIAL_FORMAT);
  vdc.Load();
  dc.SetOwnData(false);
  dc.SetMesh(vdc.GetMesh());
  for(auto it : vdc.GetFieldMap())
  {
    if(it.first.find("vol_frac_") != std::string::npos)
      dc.RegisterField(it.first, it.second);
  } 
}

void
replacementRuleTest(const std::string &shapeFile, const std::string &policyName,
  int policy, double tolerance)
{
  // Make potential baseline filenames for this test. Make a policy-specific
  // baseline that we can check first. If it is not present, the next baseline
  // is tried.
  std::string baselineName(yaml_root(shapeFile));
  std::vector<std::string> baselinePaths;
  // Example /path/to/axom/src/quest/tests/baseline/quest_intersection_shaper/cuda
  baselinePaths.push_back(pjoin(baselineDir(), policyName));
  // Example: /path/to/axom/src/quest/tests/baseline/quest_intersection_shaper/
  baselinePaths.push_back(baselineDir());

#if 0
  // Get the test info so we can use it to help construct baseline names.
  const testing::TestInfo* const test_info =
    testing::UnitTest::GetInstance()->current_test_info();
  printf("We are in test %s of test suite %s.\n",
         test_info->name(),
#endif

  // Need to make a target mesh
  SLIC_INFO(axom::fmt::format("Creating dc {}", baselineName));
  sidre::MFEMSidreDataCollection dc(baselineName, nullptr, true);
  makeTestMesh(dc);

  // Set up shapes.
  SLIC_INFO(axom::fmt::format("Reading shape set from {}", shapeFile));
  klee::ShapeSet shapeSet(klee::readShapeSet(shapeFile));

  // Need to do the pipeline of the shaping driver.
  SLIC_INFO(axom::fmt::format("Shaping materials..."));
  const int refinementLevel = 7;
#ifdef AXOM_USE_MPI
  // This has to happen here because the shaper gets its communicator from it.
  // If we do it before the mfem mesh is added to the data collection then the
  // data collection communicator gets set to MPI_COMM_NULL, which is bad for
  // the C2C reader.
  dc.SetComm(MPI_COMM_WORLD);
#endif
  quest::IntersectionShaper shaper(shapeSet, &dc);
  shaper.setLevel(refinementLevel);
  shaper.setExecPolicy(policy);

  // Borrowed from shaping_driver.
  const klee::Dimensions shapeDim = shapeSet.getDimensions();
  for(const auto& shape : shapeSet.getShapes())
  {
    SLIC_INFO(axom::fmt::format("\tshape {} -> material {}", shape.getName(), shape.getMaterial()));

    // Load the shape from file
    shaper.loadShape(shape);
    slic::flushStreams();

    // Apply the specified geometric transforms
    shaper.applyTransforms(shape);
    slic::flushStreams();

    // Generate a spatial index over the shape
    shaper.prepareShapeQuery(shapeDim, shape);
    slic::flushStreams();

    // Query the mesh against this shape
    shaper.runShapeQuery(shape);
    slic::flushStreams();

    // Apply the replacement rules for this shape against the existing materials
    shaper.applyReplacementRules(shape);
    slic::flushStreams();

    // Finalize data structures associated with this shape and spatial index
    shaper.finalizeShapeQuery();
    slic::flushStreams();
  }

#ifdef VISUALIZE_DATASETS
  saveVisIt("", baselineName, dc);
#endif
#ifdef GENERATE_BASELINES
  for(const auto &path : baselinePaths)
  {
    SLIC_INFO(axom::fmt::format("Saving baseline to {}", path));
#if 0
    saveVisIt(path, baselineName, dc);
#else
    dc.SetPrefixPath(path);
    dc.SetCycle(0);
    dc.SetPadDigits(0);
    dc.Save();
#endif
  }
#endif

// TODO: I might want an auto compare for generating baselines so I know if I need a policy-specific baseline.

  // Need to get the MFEM mesh out and compare to expected results
  bool success = false;
  for(const auto &path : baselinePaths)
  {
    try
    {
#if 0
/**
 This approach does not work because MFEM decides to mess up the filename.

MFEM Warning: Unable to open mesh file: /usr/WS2/whitlocb/Axom/axom_data/quest/regression/quest_intersection_shaper/omp/case1_000000/mesh.000000
 ... in function: void mfem::VisItDataCollection::LoadMesh()
 ... in file: fem/datacollection.cpp:574

 */

      // Load as VisIt and turn into Sidre. This is easier than getting Sidre
      // to read the data.
      mfem::VisItDataCollection vdc(baselineName);
      vdc.SetFormat(mfem::DataCollection::SERIAL_FORMAT);
      vdc.SetPrefixPath(path);
      vdc.Load();
      sidre::MFEMSidreDataCollection baselineDC(baselineName, nullptr, false);
      loadVisIt(vdc, baselineDC);
#else

/**
This approach using Sidre fails because it thinks that the volume fraction arrays
do not conform to the blueprint. Well, that has to be Sidre's fault since it wrote
the data.

WARNING in line 723 of file /usr/WS2/whitlocb/Axom/axom/src/axom/sidre/core/MFEMSidreDataCollection.cpp]
MESSAGE=MFEMSidreDataCollection blueprint verification failed:
*/
      // Try loading the baseline.
      SLIC_INFO(axom::fmt::format("Load baseline {} from {}", baselineName, path));
      sidre::MFEMSidreDataCollection baselineDC(baselineName, nullptr, true);
      baselineDC.SetPrefixPath(path);
#ifdef AXOM_USE_MPI
      baselineDC.SetComm(MPI_COMM_WORLD);
#endif
      // Make sure that the data collection does not include cycle in the name.
      baselineDC.SetPadDigits(0);
      baselineDC.Load();
#endif
      // Compare the baseline to the current DC.
      SLIC_INFO(axom::fmt::format("Comparing to baseline ", pjoin(path,baselineName)));
      success = compareDC(dc, baselineDC, tolerance);
      break;
    }
    catch(...)
    {
      SLIC_INFO(axom::fmt::format("Could not load {} from {}!", baselineName, path));
    }
  }
  EXPECT_EQ(success, true);
}

void
replacementRuleTestSet(const std::vector<std::string> &cases,
  const std::string &policyName, int policy, double tolerance)
{
  for(const auto &c : cases)
  {
    replacementRuleTest(testData(c), policyName, policy, tolerance);
  }
}

// Define testing functions for different modes.
#if defined(AXOM_SEQ_INTERSECTION_SHAPER_WORKS)
TEST(IntersectionShaperTest, case1_seq)
{
  constexpr double tolerance = 1.e-8;
  replacementRuleTestSet(case1, "seq", quest::IntersectionShaper::seq, tolerance);
}
#endif
#if defined(AXOM_USE_OPENMP)
TEST(IntersectionShaperTest, case1_omp)
{
  constexpr double tolerance = 1.e-8;
  replacementRuleTestSet(case1, "omp", quest::IntersectionShaper::omp, tolerance);
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(IntersectionShaperTest, case1_cuda)
{
  constexpr double tolerance = 1.e-8;
  replacementRuleTestSet(case1, "cuda", quest::IntersectionShaper::cuda, tolerance);
}
#endif
#if defined(AXOM_USE_HIP)
TEST(IntersectionShaperTest, case1_hip)
{
  constexpr double tolerance = 1.e-8;
  replacementRuleTestSet(case1, "hip", quest::IntersectionShaper::hip, tolerance);
}
#endif

int main(int argc, char* argv[])
{
  int result = 0;
#ifdef AXOM_USE_MPI
  // This is needed because of Axom's c2c reader.
  MPI_Init(&argc, &argv);
  // See if this aborts right away.
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#endif

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  result = RUN_ALL_TESTS();
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif
  return result;
}
