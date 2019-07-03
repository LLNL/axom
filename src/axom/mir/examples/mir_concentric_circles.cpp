// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"  // for axom macros
#include "axom/slic.hpp"
#include "axom/mir.hpp"  // for Mir classes & functions

#include <string>

// namespace aliases
namespace mir = axom::mir;

//--------------------------------------------------------------------------------

std::string usageString()
{
   return "Args are <grid size> <number of circles> <output file path>";
}

int main( int argc, char** argv )
{
   axom::slic::UnitTestLogger logger;  // create & initialize test logger
   axom::slic::setLoggingMsgLevel( axom::slic::message::Info );
  
  if (argc != 4)
  {
    SLIC_WARNING("Incorrect number of args. " << usageString() );
    return 1;
  }

  try
  {
    // Parse the command line arguments
    int gridSize = std::stoi(argv[1]);
    int numCircles = std::stoi(argv[2]);
    std::string outputFilePath = std::string(argv[3]);

    // Initialize a mesh for testing MIR
    auto timer = axom::utilities::Timer(true);
    mir::MeshTester tester;
    mir::MIRMesh testMesh = tester.initTestCaseFive(gridSize, numCircles);
    timer.stop();
    SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

    // Begin material interface reconstruction
    timer.start();
    mir::InterfaceReconstructor reconstructor;
    mir::MIRMesh processedMesh;
    reconstructor.computeReconstructedInterface(testMesh, processedMesh); 
    timer.stop();
    SLIC_INFO("Material interface reconstruction time: "
          << timer.elapsedTimeInMilliSec() << " ms.");

    // Output results
    processedMesh.writeMeshToFile(outputFilePath + "outputConcentricCircles.vtk");
    
    return 0;
  }
  catch (std::invalid_argument const &e)
  {
    SLIC_WARNING("Bad input. " << usageString() );
    return 1;
  }
  catch (std::out_of_range const &e)
  {
    SLIC_WARNING("Integer overflow. " << usageString() );
    return 1;
  }
}
