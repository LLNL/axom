/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include "gtest/gtest.h"

#include "quest/MeshTester.hpp"

#include "slic/slic.hpp"
#include "mint/Mesh.hpp"
#include "quest/STLReader.hpp"
#include "axom/path_config.hpp"
#include "axom_utils/FileUtilities.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>

std::ostream & operator<<(std::ostream & os, const std::vector<int> & m)
{
  for (unsigned int i = 0; i < m.size(); ++i) {
    os << m[i] << "  ";
  }
  return os;
}

std::ostream & operator<<(std::ostream & os,
                          const std::vector< std::pair<int, int> > & m)
{
  for (unsigned int i = 0; i < m.size(); ++i) {
    os << "(" << m[i].first << " " << m[i].second << ")  ";
  }
  return os;
}

template<typename T>
void reportVectorMismatch(const std::vector<T> & standard,
                          const std::vector<T> & result,
                          const std::string & label)
{
  std::vector<T> missing, unexpected;

  std::set_difference(standard.begin(), standard.end(),
                      result.begin(),   result.end(),
                      std::inserter(missing, missing.begin()));
  std::set_difference(result.begin(),   result.end(),
                      standard.begin(), standard.end(),
                      std::inserter(unexpected, unexpected.begin()));

  EXPECT_TRUE(missing.size() == 0) << "Missing " << missing.size() <<
    " " << label << ":" << std::endl << missing;
  EXPECT_TRUE(unexpected.size() == 0) << "Unexpectedly, " << unexpected.size() <<
    " extra " << label << ":" << std::endl << unexpected;
}

void runIntersectTest(const std::string &test,
		      const std::string &tfname,
		      const std::string &tname,
		      const std::vector< std::pair<int, int> > & expisect,
		      const std::vector< int > & expdegen)
{
  if (! axom::utilities::filesystem::pathExists(test)) {
    SLIC_INFO("Test file does not exist; reporting success and skipping: " << test);
    SUCCEED();
    return;
  }

  SCOPED_TRACE(tname);

  typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

  // read in the test file into a Mesh
  axom::quest::STLReader reader;
  reader.setFileName( tfname );
  reader.read();

  // Get surface mesh
  axom::mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader.getMesh( static_cast<TriangleMesh*>( surface_mesh ) );
  // call findTriMeshIntersections() and compare results

  std::vector< int > degenerate;
  std::vector< std::pair<int, int> > collisions = 
    axom::quest::findTriMeshIntersections(surface_mesh, degenerate);

  // report discrepancies
  // assume that expisect and expdegen are sorted (as required by
  // std::set_difference)
  std::sort(collisions.begin(), collisions.end());
  std::sort(degenerate.begin(), degenerate.end());

  reportVectorMismatch(expisect, collisions, "triangle collisions");
  reportVectorMismatch(expdegen, degenerate, "degenerate triangles");

  delete surface_mesh;
}

void splitStringToIntPairs(std::string & pairs, std::vector< std::pair<int, int> > & dat)
{
  std::istringstream iss(pairs);
  std::string a, b;
  while (std::getline(iss, a, ' ') && std::getline(iss, b, ' ')) {
    dat.push_back(std::make_pair(std::atoi(a.c_str()), std::atoi(b.c_str())));
  }
}

void splitStringToInts(std::string & pairs, std::vector< int > & dat)
{
  std::istringstream iss(pairs);
  std::string a;
  while (std::getline(iss, a, ' ')) {
    dat.push_back(std::atoi(a.c_str()));
  }
}

std::string readIntersectTest(std::string & test,
                              std::string & tfname,
                              std::vector< std::pair<int, int> > & expisect,
                              std::vector< int > & expdegen)
{
  // given a test file path in argument test,
  // return the display name for the test (from the first line of the file).
  // Output argument tfname supplies the mesh file to read in (second line, path relative to test)
  // Output arg expisect (third line) supplies the expected intersecting triangles
  // Output arg expdegen (fourth line) supplies the expected degenerate triangles

  std::string testdir;
  axom::utilities::filesystem::getDirName(testdir, test);

  std::ifstream testfile(test.c_str());
  std::string retval;
  std::getline(testfile, retval);
  std::getline(testfile, tfname);
  tfname = axom::utilities::filesystem::joinPath(testdir, tfname);
  std::string splitline;
  std::getline(testfile, splitline);
  splitStringToIntPairs(splitline, expisect);
  std::sort(expisect.begin(), expisect.end());
  std::getline(testfile, splitline);
  splitStringToInts(splitline, expdegen);
  std::sort(expdegen.begin(), expdegen.end());

  return retval;
}
                              

std::vector<std::string> findIntersectTests()
{
  std::vector<std::string> tests;

  std::string catalogue =
    axom::utilities::filesystem::joinPath(AXOM_SRC_DIR,
                                          "components/quest/data/meshtester/catalogue.txt");

  std::string testdir;
  axom::utilities::filesystem::getDirName(testdir, catalogue);

  // open file, and put each of its lines into return value tests.
  std::ifstream catfile(catalogue.c_str());
  std::string line;
  while (std::getline(catfile, line)) {
    tests.push_back(axom::utilities::filesystem::joinPath(testdir, line));
  }
  return tests;
}

TEST( quest_mesh_tester, surfacemesh_self_intersection )
{
  std::vector<std::string> tests = findIntersectTests();

  if (tests.size() < 1) {
    SLIC_INFO("*** No surface mesh self intersection tests found.");

    SUCCEED();
  }

  std::vector<std::string>::iterator it = tests.begin();
  for ( ; it != tests.end(); ++it) {
    std::string & test = *it;
    std::vector< std::pair<int, int> > expisect;
    std::vector< int > expdegen;
    std::string tfname;
    std::string tname = readIntersectTest(test, tfname, expisect, expdegen);

    // std::cout << "Read test." << std::endl <<
    //   "  Name:  " << test << std::endl <<
    //   "  File:  " << tfname << std::endl <<
    //   "  Title: " << tname << std::endl <<
    //   "  Intersecting tri pairs:  ";
    // for (int i = 0; i < expisect.size(); ++i) {
    //   std::cout << expisect[i].first << " " << expisect[i].second << "  ";
    // }
    // std::cout << std::endl << "  Degenerate tris:  ";
    // for (int i = 0; i < expdegen.size(); ++i) {
    //   std::cout << expdegen[i] << " ";
    // }
    // std::cout << std::endl;

    runIntersectTest(test, tfname, tname, expisect, expdegen);
  }
}
