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
namespace numerics = axom::numerics;
namespace slam = axom::slam;
namespace mir = axom::mir;
namespace fs = axom::utilities::filesystem;

//--------------------------------------------------------------------------------

enum InputStatus
{
  SUCCESS,
  INVALID,
  SHOWHELP
};

struct Input
{
   int m_test_case; // valid values 1,2,3,4,5
   bool m_should_iterate;
   int m_iter_count;
   double m_iter_percent;
   bool m_verbose;
   std::string m_output_dir;
   InputStatus m_status;

   Input()
      : m_test_case(2)
      , m_should_iterate(false)
      , m_iter_count(100)
      , m_iter_percent(.3)
      , m_verbose(false)
      , m_status(SUCCESS)
   {
      m_output_dir = fs::joinPath(AXOM_BIN_DIR, "mir_examples");
   }

   Input(int argc, char** argv) : Input()
   {
      for (int i = 1 ; i < argc ; /* increment i in loop */)
      {
        std::string arg = argv[i];
        if (arg == "--test-case")
        {
          m_test_case = std::stoi(argv[++i]);
        }
        else if (arg == "--output-dir")
        {
           m_output_dir = argv[++i];
        }
        else if (arg == "--iter-count")
        {
           m_iter_count = std::stoi(argv[++i]);
           m_should_iterate = true;
        }
        else if (arg == "--iter-percent")
        {
           m_iter_percent = std::stod(argv[++i]);
           m_should_iterate = true;
        }
        else if (arg == "--verbose")
        {
           m_verbose = true;
        }
        else // help or unknown parameter
        {
          if(arg != "--help" && arg != "-h")
          {
            SLIC_WARNING("Unrecognized parameter: " << arg);
            m_status = INVALID;
          }
          else
          {
             m_status = SHOWHELP;
          }
          return;
        }
        ++i;
      }

      checkTestCase();
      checkOutputDir();
      checkIterationParams();
   }


   bool shouldIterate() const { return m_should_iterate; }
   int numIterations() const { return m_iter_count; }
   int iterPercentage() const { return m_iter_percent; }


   void showhelp()
   {
     std::cout
       << "Argument usage:"
        "\n  --help            Show this help message."
        "\n  --test-case N     Mesh test case.  Default N = 2."
        "\n                    Valid values {1,2,3,4,5,6}"
        "\n  --output-dir dir  Directory for output mesh"
        "\n                    Default is: '${AXOM_BIN_DIR}/mir_examples'"
        "\n  --iter-count N    Number of iterations for iterative algorithm"
        "\n                    Defaults to 100."
        "\n                    Setting a value triggers iterative algorithm"
        "\n  --iter-percent D  Volume diff percentage for iterative algorithm"
        "\n                    Must be between 0 and 1. Defaults to 0.3"
        "\n                    Setting a value triggers iterative algorithm"
        "\n  --verbose         Increases verbosity of output"
       << std::endl << std::endl;
   };

private:
   void checkTestCase()
   {
      if(m_test_case < 1 || m_test_case > 6)
      {
         m_status = INVALID;
         SLIC_WARNING("Invalid test case " << m_test_case);
      }
   }

   void checkOutputDir()
   {
      if(! fs::pathExists(m_output_dir))
      {
         fs::makeDirsForPath(m_output_dir);
      }
   }

   void checkIterationParams()
   {
      if(m_should_iterate)
      {
         if(m_iter_count < 1)
         {
            m_status = INVALID;
            SLIC_WARNING("Invalid iteration count " << m_iter_count);
         }

         if(m_iter_percent <= 0. || m_iter_percent > 1.)
         {
            m_status = INVALID;
            SLIC_WARNING("Invalid iteration percentage " << m_iter_percent);
         }
      }
   }

};

/*!
 * \brief Tutorial main showing how to initialize test cases and perform mir.
 */
int main(int argc, char** argv)
{
  axom::slic::UnitTestLogger logger;  // create & initialize test logger
  axom::slic::setLoggingMsgLevel( axom::slic::message::Info );
  
  // Parse arguments
  Input params(argc, argv);

  if (params.m_status != SUCCESS)
  {
    if (params.m_status == SHOWHELP)
    {
      params.showhelp();
      return 0;
    }
    else if (params.m_status == INVALID)
    {
      params.showhelp();
      return 1;
    }
  }

  mir::MIRMesh testMesh;
  mir::MeshTester tester;

  auto timer = axom::utilities::Timer(true);

  switch(params.m_test_case)
  {
  case 1: testMesh = tester.initTestCaseOne(); break;
  case 2: testMesh = tester.initTestCaseTwo(); break;
  case 3: testMesh = tester.initTestCaseThree(); break;
  case 4: testMesh = tester.initTestCaseFour(); break;
  case 5: testMesh = tester.initTestCaseFive(25, 12); break;
  case 6: testMesh = tester.initTestCaseSix(15, 3); break;
  }

  timer.stop();
  SLIC_INFO( "Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

  SLIC_INFO("Test mesh is " << (testMesh.isValid(true) ? "" : " NOT") << " valid.");

  if(params.m_verbose)
  {
     SLIC_INFO("Initial mesh:");
     testMesh.print();
  }

  // Begin material interface reconstruction
  timer.start();

  mir::MIRMesh processedMesh;
  mir::InterfaceReconstructor reconstructor;
  
  if(! params.shouldIterate())  // Process once, with original Meredith algorithm
  {
     reconstructor.computeReconstructedInterface(testMesh, processedMesh);
  }
  else // use iterative algorithm
  {
     int n = params.numIterations();
     double p = params.iterPercentage();

     reconstructor.computeReconstructedInterfaceIterative(testMesh, n, p, processedMesh);
  }

  timer.stop();
  SLIC_INFO( "Reconstruction time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Output results
  processedMesh.writeMeshToFile(params.m_output_dir, "processedMesh.vtk");

  using VolFracs =  std::vector<std::vector<axom::float64> >;
  timer.start();
  VolFracs materialVolumeFractionsElement = processedMesh.computeOriginalElementVolumeFractions();
  timer.stop();
  SLIC_INFO( "Computing volumes took: " << timer.elapsedTimeInMilliSec() << " ms.");

  if(params.m_verbose)
  {
     SLIC_INFO("Final mesh:");
     processedMesh.print();
  }

  return 0;
}

//--------------------------------------------------------------------------------
