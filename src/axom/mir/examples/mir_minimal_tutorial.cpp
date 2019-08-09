// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

// _mir_header_start
#include "axom/mir.hpp"
// _mir_header_end

#include <string>

// namespace aliases
namespace mir = axom::mir;
namespace fs = axom::utilities::filesystem;

/*!
 * \brief Tutorial main showing how to initialize test cases and perform mir.
 */
int main(int argc, char** argv)
{
  // _mir_main_loop_start
  // Initialize the mesh
  mir::MIRMesh testMesh;
  mir::MeshTester tester;
  testMesh = tester.initTestCaseOne();

  // Perform material interface reconstruction
  mir::MIRMesh processedMesh;
  mir::InterfaceReconstructor reconstructor;
  reconstructor.computeReconstructedInterface( testMesh, processedMesh );

  // Output the results to a .vtk file
  std::string outputDirectory = fs::joinPath(AXOM_BIN_DIR, "mir_examples/processedMesh.vtk");
  processedMesh.writeMeshToFile( outputDirectory );
  // _mir_main_loop_end
  return 0;
}