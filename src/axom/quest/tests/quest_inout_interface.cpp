// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"  // for gtest macros

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "axom/quest/interface/inout.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"  // for test that reads
                                                           // in a mesh without
                                                           // the interface
#include <string>

/// Test fixture for quest::inout_query interface
class InOutInterfaceTest : public ::testing::Test
{
public:
  static const int DIM = 3;

  using CoordType = double;
  using InOutPoint = axom::primal::Point<CoordType, DIM>;
  using InOutBBox = axom::primal::BoundingBox<CoordType, DIM>;

protected:
  virtual void SetUp()
  {
    // Setup meshfile
#ifdef AXOM_DATA_DIR
    namespace fs = axom::utilities::filesystem;
    const std::string DATA_DIR = fs::joinPath(AXOM_DATA_DIR, "quest");
    const std::string fileName = "sphere.stl";
    meshfile = fs::joinPath(DATA_DIR, fileName);
#else
    FAIL() << "quest_inout_interface test requires AXOM_DATA_DIR";
#endif
  }

  std::string meshfile;
};

TEST_F(InOutInterfaceTest, initialize_and_finalize)
{
  EXPECT_TRUE(axom::utilities::filesystem::pathExists(this->meshfile));

  // InOut begins uninitialized
  EXPECT_FALSE(axom::quest::inout_initialized());

  // Initialize the InOut query
  EXPECT_EQ(0, axom::quest::inout_init(this->meshfile));

  // InOut should now be initialized
  EXPECT_TRUE(axom::quest::inout_initialized());

  // Finalize the InOut query
  EXPECT_EQ(0, axom::quest::inout_finalize());

  // InOut should no longer  be initialized
  EXPECT_FALSE(axom::quest::inout_initialized());
}

TEST_F(InOutInterfaceTest, logger_inited)
{
  const bool origSlicInited = axom::slic::isInitialized();

  auto origLogLevel = axom::slic::getLoggingMsgLevel();

  // Initialize the InOut query
  EXPECT_EQ(0, axom::quest::inout_init(this->meshfile));

  // slic is initialized now
  EXPECT_TRUE(axom::slic::isInitialized());

  // Finalize the InOut query
  EXPECT_EQ(0, axom::quest::inout_finalize());

  // After finalizing, slic should be in the same status that it was originally
  EXPECT_EQ(origSlicInited, axom::slic::isInitialized());

  if(origSlicInited)
  {
    EXPECT_EQ(origLogLevel, axom::slic::getLoggingMsgLevel());
  }
}

TEST_F(InOutInterfaceTest, initialize_from_mesh)
{
  EXPECT_TRUE(axom::utilities::filesystem::pathExists(this->meshfile));

  axom::mint::Mesh* mesh = nullptr;
  int rc = axom::quest::internal::read_mesh(this->meshfile, mesh);
  EXPECT_EQ(0, rc);

  // Initialize the InOut query
  EXPECT_EQ(0, axom::quest::inout_init(mesh));

  // InOut should now be initialized
  EXPECT_TRUE(axom::quest::inout_initialized());

  // Finalize the InOut query
  EXPECT_EQ(0, axom::quest::inout_finalize());

  // InOut should no longer  be initialized
  EXPECT_FALSE(axom::quest::inout_initialized());
}

TEST_F(InOutInterfaceTest, query_properties)
{
  const int failCode = axom::quest::QUEST_INOUT_FAILED;
  const int successCode = axom::quest::QUEST_INOUT_SUCCESS;

  using PointType = InOutInterfaceTest::InOutPoint;
  using BBoxType = InOutInterfaceTest::InOutBBox;
  PointType lo, hi, cm;

  // first, test before initializing
  {
    EXPECT_FALSE(axom::quest::inout_initialized());

    // The following should return an error since initialized is false
    SLIC_INFO("--[==[ \n"
              << "\t The following three calls might emit warning messages."
              << " This is expected.\n");

    EXPECT_EQ(failCode, axom::quest::inout_mesh_min_bounds(lo.data()));
    EXPECT_EQ(failCode, axom::quest::inout_mesh_max_bounds(hi.data()));
    EXPECT_EQ(failCode, axom::quest::inout_mesh_center_of_mass(cm.data()));
    EXPECT_EQ(failCode, axom::quest::inout_get_dimension());
    SLIC_INFO("--]==]");
  }

  // next, test after initialization
  {
    axom::quest::inout_init(this->meshfile);

    EXPECT_EQ(successCode, axom::quest::inout_mesh_min_bounds(lo.data()));
    EXPECT_EQ(successCode, axom::quest::inout_mesh_max_bounds(hi.data()));
    EXPECT_EQ(successCode, axom::quest::inout_mesh_center_of_mass(cm.data()));

    const int expSpatialDim = 3;
    int spatialDim = axom::quest::inout_get_dimension();
    EXPECT_EQ(expSpatialDim, spatialDim);

    SLIC_INFO("Mesh bounding box is " << BBoxType(lo, hi));
    SLIC_INFO("Mesh center of mass is " << cm);
    SLIC_INFO("Spatial dimension of query is " << spatialDim);

    axom::quest::inout_finalize();
  }
}

TEST_F(InOutInterfaceTest, set_params)
{
  const int failCode = axom::quest::QUEST_INOUT_FAILED;
  const int successCode = axom::quest::QUEST_INOUT_SUCCESS;
  const double EPS = 1E-6;

  // Set parameters at the appropriate time
  {
    EXPECT_FALSE(axom::quest::inout_initialized());

    EXPECT_EQ(successCode, axom::quest::inout_set_verbose(true));
    EXPECT_EQ(successCode, axom::quest::inout_set_vertex_weld_threshold(EPS));
  }

  // Initialize the query
  {
    SLIC_INFO("*** About to initialize with verbose output.");
    axom::quest::inout_init(this->meshfile);
    SLIC_INFO("*** End verbose output.\n");
  }

  // Set parameters at an inappropriate time
  {
    SLIC_INFO("*** The following calls might emit warning messages.");
    SLIC_INFO("--[==[");

    EXPECT_EQ(failCode, axom::quest::inout_set_verbose(true));
    EXPECT_EQ(failCode, axom::quest::inout_set_vertex_weld_threshold(EPS));

    SLIC_INFO("--]==]");
  }

  axom::quest::inout_finalize();
}

TEST_F(InOutInterfaceTest, query)
{
  typedef InOutInterfaceTest::InOutPoint PointType;
  PointType query;

  axom::quest::inout_init(this->meshfile);

  // test an inside point
  EXPECT_TRUE(axom::quest::inout_evaluate(0, 0, 0));

  // test an outside point
  EXPECT_FALSE(axom::quest::inout_evaluate(10, 10, 10));

  axom::quest::inout_finalize();
}

int main(int argc, char** argv)
{
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
#endif
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  int result = RUN_ALL_TESTS();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return result;
}
