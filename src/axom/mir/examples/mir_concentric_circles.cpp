// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"  // for axom macros
#include "axom/mir.hpp"  // for Mir classes & functions
#include "axom/slam.hpp"

#include <chrono>
#include <string>

using Clock = std::chrono::high_resolution_clock;

// namespace aliases
namespace numerics = axom::numerics;
namespace slam = axom::slam;
namespace mir = axom::mir;

//--------------------------------------------------------------------------------

/*!
 * \brief Tutorial main
 */
int main( int argc, char** argv )
{
  
  if (argc != 4)
  {
    printf("Incorrect number of args. Args are <grid size> <number of circles> <output file path>\n");
    return 0;
  }

  try
  {
    // Parse the command line arguments
    int gridSize = std::stoi(argv[1]);
    int numCircles = std::stoi(argv[2]);
    std::string outputFilePath = std::string(argv[3]);

    // Intialize a mesh for testing MIR
    auto startTime = Clock::now();
    mir::MeshTester tester;
    mir::MIRMesh testMesh = tester.initTestCaseFive(gridSize, numCircles);
    auto endTime = Clock::now();
    std::cout << "Mesh init time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;

    // Begin material interface reconstruction
    startTime = Clock::now();
    mir::InterfaceReconstructor reconstructor;
    mir::MIRMesh processedMesh = reconstructor.computeReconstructedInterface(testMesh); 
    endTime = Clock::now();
    std::cout << "Material interface reconstruction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;

    // Output results
    processedMesh.writeMeshToFile(outputFilePath + "outputConcentricCircles.vtk");
    
    return 0;
  }
  catch (std::invalid_argument const &e)
  {
    printf("Bad input. Arguments are <grid size> <number of circles> <output file path>\n");
    return 0;
  }
  catch (std::out_of_range const &e)
  {
    printf("Integer overflow. Arguments are <grid size> <number of circles> <output file path>\n");
    return 0;
  }
}