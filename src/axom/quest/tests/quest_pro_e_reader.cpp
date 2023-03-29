// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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
  const std::string INVALID_FILE = "foo.creo";
  axom::quest::ProEReader reader;
  reader.setFileName(INVALID_FILE);
  int status = reader.read();
  EXPECT_TRUE(status != 0);
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader, read_to_invalid_mesh)
{
  const char* IGNORE_OUTPUT = ".*";
  const std::string filename = "tet.creo";

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

  const std::string filename = "tet.creo";

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
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}