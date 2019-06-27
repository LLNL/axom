// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"  // for axom macros
#include "axom/mir.hpp"  // for Mir classes & functions
#include "axom/slam.hpp"

#include <chrono>
using Clock = std::chrono::high_resolution_clock;

// namespace aliases
namespace numerics = axom::numerics;
namespace slam = axom::slam;
namespace mir = axom::mir;

//--------------------------------------------------------------------------------

/*!
 * \brief Tutorial main showing how to initialize test cases and perform mir.
 */
int main( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
{
  
  // Intialize a mesh for testing MIR
  auto startTime = Clock::now();
  mir::MeshTester tester;
  // mir::MIRMesh testMesh = tester.initTestCaseOne();
  mir::MIRMesh testMesh = tester.initTestCaseTwo();
  // mir::MIRMesh testMesh = tester.initTestCaseThree();
  // mir::MIRMesh testMesh = tester.initTestCaseFour();
  // mir::MIRMesh testMesh = tester.initTestCaseFive(25, 12);

  auto endTime = Clock::now();
  std::cout << "Mesh init time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;

  // Begin material interface reconstruction
  startTime = Clock::now();

  mir::MIRMesh processedMesh;
  
  mir::InterfaceReconstructor reconstructor;
  reconstructor.computeReconstructedInterface(testMesh, processedMesh);                // Process once, with original Meredith algorithm
  // reconstructor.computeReconstructedInterfaceIterative(testMesh, 100, 0.3, processedMesh);  // 100 iterations, 20 percent with iterative Meredith algorithm

  endTime = Clock::now();
  std::cout << "Material interface reconstruction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;

  // Output results
  processedMesh.writeMeshToFile("/Users/sterbentz3/Desktop/testOutputNewMARKONEW4.vtk");

  std::vector<std::vector<axom::float64> > materialVolumeFractionsElement = processedMesh.computeOriginalElementVolumeFractions();
  
    // // Print out the results 
    // printf("elementVolumeFractions: {\n");
    // for (int matID = 0; matID < processedMesh.m_numMaterials; ++matID)
    // {
    //   printf("Material %d: {", matID);
    //   for (int eID = 0; eID < materialVolumeFractionsElement[matID].size(); ++eID)
    //   {
    //     printf(" %f,", materialVolumeFractionsElement[matID][eID]);
    //   }
    //   printf("}\n");
    // }
    // printf("}\n");
    // // END TESTING

  return 0;
}

//--------------------------------------------------------------------------------