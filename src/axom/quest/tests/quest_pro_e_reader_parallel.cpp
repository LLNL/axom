// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/slic.hpp"

#include "axom/quest/readers/PProEReader.hpp"

#include "gtest/gtest.h"

#include "mpi.h"

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
TEST(quest_pro_e_reader_parallel, missing_file)
{
  const std::string INVALID_FILE = "foo.creo";
  axom::quest::PProEReader reader(MPI_COMM_WORLD);
  reader.setFileName(INVALID_FILE);
  int status = reader.read();
  EXPECT_TRUE(status == -1);
}

//------------------------------------------------------------------------------
TEST(quest_pro_e_reader_parallel, read_file)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const double x_expected[] = {-1.0, 1.0, 0.0, 0.0};
  const double y_expected[] = {0.0, 0.0, 1.0, 0.0};
  const double z_expected[] = {0.0, 0.0, 0.0, 1.0};

  const std::string filename = "tet.creo";

  // STEP 0: generate a temporary Pro/E file for testing
  if(rank == 0)
  {
    generate_pro_e_file(filename);
  }

  // STEP 1: create an Pro/E reader and read-in the mesh data
  axom::quest::PProEReader reader(MPI_COMM_WORLD);
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

  // // STEP 4: remove temporary Pro/E file
  if(rank == 0)
  {
    std::remove(filename.c_str());
  }
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  MPI_Init(&argc, &argv);

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}