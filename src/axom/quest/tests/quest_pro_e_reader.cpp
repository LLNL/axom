// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/mint/utils/vtk_utils.hpp"  // for write_vtk
#include "axom/quest/readers/ProEReader.hpp"
#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

#include <fstream>

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
/*!
 * \brief Generates a Pro/E file consisting of a single tetrahedron
 * \param [in] file the name of the file to generate.
 * \pre file.empty() == false
 */
void generate_pro_e_file(const std::string& file)
{
  EXPECT_FALSE(file.empty());

  std::ofstream ofs(file.c_str());
  EXPECT_TRUE(ofs.is_open());

  ofs << "# Comment header to ignore" << std::endl;
  ofs << "# Another comment" << std::endl;
  ofs << "4 1" << std::endl;
  ofs << "node1 -1.0 0.0 0.0" << std::endl;
  ofs << "node2 1.0 0.0 0.0" << std::endl;
  ofs << "node3 0.0 1.0 0.0" << std::endl;
  ofs << "node4 0.0 0.0 1.0" << std::endl;
  ofs << "tet1 node1 node2 node3 node4" << std::endl;

  ofs.close();
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_missing_file)
{
  const std::string INVALID_FILE = "foo.proe";
  axom::quest::ProEReader reader;
  reader.setFileName(INVALID_FILE);
  int status = reader.read();
  EXPECT_TRUE(status != 0);
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_to_invalid_mesh)
{
  const char* IGNORE_OUTPUT = ".*";
  const std::string filename = "tet.proe";

  // STEP 0: generate a temporary Pro/E file for testing
  generate_pro_e_file(filename);

  // STEP 1: constructs mesh object to read in the mesh to
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> trimesh(
    2,
    axom::mint::TRIANGLE);
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> hexmesh(3,
                                                                 axom::mint::HEX);

  // STEP 2: read in the Pro/E mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(filename);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 3: death tests

  // read the Pro/E mesh data to a 2D axom::mint::Mesh should fail
  EXPECT_DEATH_IF_SUPPORTED(reader.getMesh(&trimesh), IGNORE_OUTPUT);

  // read the Pro/E mesh data to a axom::mint::Mesh that has a different cell type
  EXPECT_DEATH_IF_SUPPORTED(reader.getMesh(&hexmesh), IGNORE_OUTPUT);

  // STEP 4: remove Pro/E file
  std::remove(filename.c_str());
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e)
{
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  const std::string filename = "tet.proe";

  // STEP 0: generate a temporary Pro/E file for testing
  generate_pro_e_file(filename);

  // STEP 1: create an Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(filename);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), 1);
  EXPECT_EQ(mesh.getNumberOfNodes(), 4);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode],
                x_expected[inode],
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode],
                y_expected[inode],
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode],
                z_expected[inode],
                std::numeric_limits<double>::epsilon());
  }  // END for all nodes

  // STEP 4: remove temporary Pro?E file
  std::remove(filename.c_str());
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_pro_e_external)
{
  constexpr axom::IndexType N_NODES = 4;
  constexpr axom::IndexType N_TETS = 1;
  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  double xin[] = {-1.0, -1.0, -1.0, -1.0};
  double yin[] = {-1.0, -1.0, -1.0, -1.0};
  double zin[] = {-1.0, -1.0, -1.0, -1.0};

  axom::IndexType conn[] = {-1, -1, -1, -1};

  const std::string filename = "tet.proe";

  // STEP 0: generate a temporary Pro/E file for testing
  generate_pro_e_file(filename);

  // STEP 1: create a Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(filename);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>
    mesh(axom::mint::TET, N_TETS, conn, N_NODES, xin, yin, zin);
  EXPECT_EQ(mesh.getNumberOfCells(), N_TETS);
  EXPECT_EQ(mesh.getNumberOfNodes(), N_NODES);

  reader.getMesh(&mesh);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  axom::IndexType numNodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    EXPECT_NEAR(x[inode],
                x_expected[inode],
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(y[inode],
                y_expected[inode],
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(z[inode],
                z_expected[inode],
                std::numeric_limits<double>::epsilon());
  }  // END for all nodes

  // STEP 4: remove temporary Pro/E file
  std::remove(filename.c_str());
}

//------------------------------------------------------------------------------
#ifdef AXOM_DATA_DIR
TEST(quest_pro_e_reader, cup_pro_e)
{
  constexpr int NUM_TETS = 574;
  constexpr double EPS = std::numeric_limits<double>::epsilon();

  // STEP 0: Get Pro/E cup example file for testing
  namespace fs = axom::utilities::filesystem;
  std::string cup = fs::joinPath(AXOM_DATA_DIR, "quest/cup.proe");

  // STEP 1: create a Pro/E reader and read-in the mesh data
  axom::quest::ProEReader reader;
  reader.setFileName(cup);
  int status = reader.read();
  EXPECT_EQ(status, 0);

  // STEP 2: reading the Pro/E mesh data into a axom::mint::Mesh
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(3, axom::mint::TET);
  reader.getMesh(&mesh);

  // STEP 3: ensure the mesh is what is expected
  EXPECT_EQ(mesh.getNumberOfCells(), NUM_TETS);
  EXPECT_EQ(mesh.getNumberOfNodes(), NUM_TETS * 4);

  const double* x = mesh.getCoordinateArray(axom::mint::X_COORDINATE);
  const double* y = mesh.getCoordinateArray(axom::mint::Y_COORDINATE);
  const double* z = mesh.getCoordinateArray(axom::mint::Z_COORDINATE);
  EXPECT_TRUE(x != nullptr);
  EXPECT_TRUE(y != nullptr);
  EXPECT_TRUE(z != nullptr);

  // STEP 4: Verify a few tetrahedra are as expected

  // Check first tetrahedron
  EXPECT_NEAR(x[0], -71.842646, EPS);
  EXPECT_NEAR(x[1], 18.091620, EPS);
  EXPECT_NEAR(x[2], 18.091620, EPS);
  EXPECT_NEAR(x[3], 16.637112, EPS);

  EXPECT_NEAR(y[0], -142.432949, EPS);
  EXPECT_NEAR(y[1], -100.389772, EPS);
  EXPECT_NEAR(y[2], -158.496833, EPS);
  EXPECT_NEAR(y[3], -116.282419, EPS);

  EXPECT_NEAR(z[0], 0.000000, EPS);
  EXPECT_NEAR(z[1], 0.000000, EPS);
  EXPECT_NEAR(z[2], 0.000000, EPS);
  EXPECT_NEAR(z[3], -53.197840, EPS);

  // Check last tetrahedron
  constexpr int LAST_OFFSET = (NUM_TETS - 1) * 4;

  EXPECT_NEAR(x[LAST_OFFSET], 128.940529, EPS);
  EXPECT_NEAR(x[LAST_OFFSET + 1], 68.261330, EPS);
  EXPECT_NEAR(x[LAST_OFFSET + 2], 121.442679, EPS);
  EXPECT_NEAR(x[LAST_OFFSET + 3], 68.137155, EPS);

  EXPECT_NEAR(y[LAST_OFFSET], -93.930323, EPS);
  EXPECT_NEAR(y[LAST_OFFSET + 1], -41.261753, EPS);
  EXPECT_NEAR(y[LAST_OFFSET + 2], -31.346653, EPS);
  EXPECT_NEAR(y[LAST_OFFSET + 3], -41.466486, EPS);

  EXPECT_NEAR(z[LAST_OFFSET], -145.985725, EPS);
  EXPECT_NEAR(z[LAST_OFFSET + 1], -122.672539, EPS);
  EXPECT_NEAR(z[LAST_OFFSET + 2], -174.327849, EPS);
  EXPECT_NEAR(z[LAST_OFFSET + 3], -173.363017, EPS);

  // Check middle tetrahedron
  constexpr int MID_OFFSET = ((NUM_TETS / 2) - 1) * 4;

  EXPECT_NEAR(x[MID_OFFSET], 15.508351, EPS);
  EXPECT_NEAR(x[MID_OFFSET + 1], -30.194340, EPS);
  EXPECT_NEAR(x[MID_OFFSET + 2], -35.211560, EPS);
  EXPECT_NEAR(x[MID_OFFSET + 3], 14.765131, EPS);

  EXPECT_NEAR(y[MID_OFFSET], 124.048744, EPS);
  EXPECT_NEAR(y[MID_OFFSET + 1], 109.246291, EPS);
  EXPECT_NEAR(y[MID_OFFSET + 2], 117.715966, EPS);
  EXPECT_NEAR(y[MID_OFFSET + 3], 78.384516, EPS);

  EXPECT_NEAR(z[MID_OFFSET], -179.442905, EPS);
  EXPECT_NEAR(z[MID_OFFSET + 1], -112.054428, EPS);
  EXPECT_NEAR(z[MID_OFFSET + 2], -178.748301, EPS);
  EXPECT_NEAR(z[MID_OFFSET + 3], -163.120569, EPS);

  // Step 5: Dump mesh file
  axom::mint::write_vtk(&mesh, "cup.vtk");
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}
