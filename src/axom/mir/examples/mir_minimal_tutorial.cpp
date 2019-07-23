// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/mir.hpp"

#include <string>

// namespace aliases
namespace mir = axom::mir;

/*!
 * \brief Tutorial main showing how to initialize test cases and perform mir.
 */
int main(int argc, char** argv)
{
  // Initialize the mesh
  mir::MIRMesh testMesh;
  mir::MeshTester tester;
  testMesh = tester.initTestCaseOne();

  // Perform material interface reconstruction
  mir::MIRMesh processedMesh;
  mir::InterfaceReconstructor reconstructor;
  reconstructor.computeReconstructedInterface(testMesh, processedMesh);

  // Output the results to a .vtk file
  processedMesh.writeMeshToFile("mir_examples", "processedMesh.vtk");

  return 0;
}