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

  // ofs << "solid triangle" << std::endl;
  // ofs << "\t facet normal 0.0 0.0 1.0" << std::endl;
  // ofs << "\t\t outer loop" << std::endl;
  // ofs << "\t\t\t vertex 0.0 0.0 0.0" << std::endl;
  // ofs << "\t\t\t vertex 1.0 0.0 0.0" << std::endl;
  // ofs << "\t\t\t vertex 0.0 1.0 0.0" << std::endl;
  // ofs << "\t\t endloop" << std::endl;
  // ofs << "\t endfacet" << std::endl;
  // ofs << "endsolid triangle" << std::endl;

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
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  return RUN_ALL_TESTS();
}