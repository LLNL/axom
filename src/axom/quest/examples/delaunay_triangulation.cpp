// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file delaunay_triangulation.cpp
 * \brief This file contains an example that creates Delaunay triangulation
 * sequentially from a set of random points contained within a rectangle
 * centered at the origin.
 * This example can run in 2D with triangles or in 3D with tetrahedrons.
 * It uses SLAM IA, a topological mesh data structure.
 */

#include "axom/slic.hpp"
#include "axom/quest.hpp"

#include <ctime>
#include <iomanip>

using namespace axom;

//This is the length of the bounding rectangle to contain the random points
double BOUNDING_BOX_DIM[] = {2.0, 2.0, 2.0};

using IndexType = int;
using DataType = double;
using Point2D = axom::primal::Point2D;
using Point3D = axom::primal::Point3D;

enum InputStatus
{
  SUCCESS,
  SHOWHELP,
  CANTOPENFILE
};

struct Input
{
  std::string inputFile;
  std::string outputVTKFile;
  int numRandPoints;
  bool outputSteps;
  int numOutputSteps;
  int dimension;

  InputStatus errorCode;

  Input()
    : outputVTKFile("delaunay_output")
    , numRandPoints(10)
    , outputSteps(true)
    , errorCode(SUCCESS) {};

  Input(int argc, char** argv);

  void showhelp()
  {
    std::cout << "Argument usage: \n"
                 "  --help(-h)             Show this help message. \n"
                 "  --nrandpt(-n) [N]      The number of random points.  "
                 "Default is N = 20.  \n"
                 "  --dim(-d) [D]          The dimension of the mesh. 2 = "
                 "triangles in 2D, 3 = tetrahedrons in 3D.\n"
                 "  --outfile(-o) [fname]  The VTK output file. Default fname "
                 "= delaunay_output.\n"
                 "  --outputsteps(-s) [S]  S number of intermediate steps will "
                 "be written out as VTK files.\n"
                 "                         S must be a strictly positive "
                 "integer, or omit to write one file per step.\n"
                 "                         By default no internediate VTK file "
                 "is written.\n"
              << std::endl;
  };
};

Input::Input(int argc, char** argv)
  : outputVTKFile("delaunay_output")
  , numRandPoints(20)
  , outputSteps(false)
  , numOutputSteps(-1)
  , dimension(2)
  , errorCode(SUCCESS)
{
  if(argc >= 2)
  {
    std::string help = argv[1];
    if(help == "--help" || help == "-h")
    {
      errorCode = SHOWHELP;
      return;
    }
    for(int i = 1; i < argc; /* increment i in loop */)
    {
      std::string arg = argv[i];
      if(arg == "--nrandpt" || arg == "-n")
      {
        numRandPoints = atoi(argv[++i]);
      }
      else if(arg == "--outputsteps" || arg == "-s")
      {
        outputSteps = true;
        if(argc > i + 1)
        {
          numOutputSteps = atoi(argv[i + 1]);
          //Check that it is actually a non-zero number
          if(numOutputSteps != 0 || *argv[i + 1] == '0')
          {
            // it's a number
            i++;
          }
          else
          {
            numOutputSteps = -1;
          }
        }
      }
      else if(arg == "--outfile" || arg == "-o")
      {
        outputVTKFile = argv[++i];
      }
      else if(arg == "--dim" || arg == "-d")
      {
        dimension = atoi(argv[++i]);
      }
      else
      {
        errorCode = SHOWHELP;
        return;
      }
      ++i;
    }
  }

  // if -s is specified but the number is obmitted, then output one file for
  // each step.
  if(outputSteps && (numOutputSteps == -1 || numOutputSteps > numRandPoints))
  {
    numOutputSteps = numRandPoints;
  }

  SLIC_INFO("Using parameter values: "
            << "\n nrandpt = " << numRandPoints << "\n dimension = " << dimension
            << "\n outputSteps = " << (numOutputSteps == -1 ? 0 : numOutputSteps)
            << "\n outfile = " << outputVTKFile << std::endl);
}

template <unsigned int dimension>
int run_delaunay(int numPoints, int numOutputVTKsteps, std::string outputVTKFile)
{
  // Create Delaunay Triangulation on random points
  using Delaunay = axom::quest::Delaunay<dimension>;
  using PointType = typename Delaunay::PointType;
  using BoundingBox = axom::primal::BoundingBox<DataType, dimension>;
  Delaunay dt;

  int stepOutputCount = 1;

  double bb_min_arr[] = {-BOUNDING_BOX_DIM[0] / 2.0,
                         -BOUNDING_BOX_DIM[1] / 2.0,
                         -BOUNDING_BOX_DIM[2] / 2.0};

  double bb_max_arr[] = {BOUNDING_BOX_DIM[0] / 2.0,
                         BOUNDING_BOX_DIM[1] / 2.0,
                         BOUNDING_BOX_DIM[2] / 2.0};

  PointType bb_min(bb_min_arr);
  PointType bb_max(bb_max_arr);
  BoundingBox bb(bb_min, bb_max);
  dt.initializeBoundary(bb);

  //Adding a number of random points
  time_t rseed_val = std::time(NULL);
  SLIC_INFO("Using " << rseed_val << " as seed for rand()");
  srand(rseed_val);

  for(int pt_i = 0; pt_i < numPoints; pt_i++)
  {
    double pt_coord[dimension];
    for(unsigned int d = 0; d < dimension; d++)
      pt_coord[d] = (double)rand() / (double)RAND_MAX * BOUNDING_BOX_DIM[d] -
        BOUNDING_BOX_DIM[d] / 2.0;
    PointType new_pt(pt_coord, dimension);

    dt.insertPoint(new_pt);

    if(numOutputVTKsteps != -1)
    {
      if(stepOutputCount * numPoints / (numOutputVTKsteps + 1) == pt_i)
      {
        std::ostringstream stm;
        stm << outputVTKFile << "_" << std::setfill('0') << std::setw(6) << pt_i
            << ".vtk";
        dt.writeToVTKFile(stm.str());
        stepOutputCount++;
      }
    }
  }

  dt.removeBoundary();  //Remove the starting rectangular box

  // Write the final mesh to a vtk file
  std::ostringstream stm;
  stm << outputVTKFile << ".vtk";
  dt.writeToVTKFile(stm.str());

  SLIC_INFO("Done!");

  return 0;
}

int main(int argc, char** argv)
{
  // Initialize the SLIC logger
  axom::slic::SimpleLogger logger;

  // Initialize default parameters and update with command line arguments:
  Input params(argc, argv);

  if(params.errorCode != SUCCESS)
  {
    if(params.errorCode == SHOWHELP)
    {
      params.showhelp();
      return EXIT_SUCCESS;
    }
    else
    {
      std::cerr << "Unknown error " << (int)params.errorCode
                << " while parsing arguments." << std::endl;
      return EXIT_FAILURE;
    }
  }

  const int numOutputSteps = params.numOutputSteps;
  const int numPoints = params.numRandPoints;
  const std::string outputVTKFile = params.outputVTKFile;
  const int dimension = params.dimension;

  if(dimension == 2)
  {
    run_delaunay<2>(numPoints, numOutputSteps, outputVTKFile);
  }
  else if(dimension == 3)
  {
    run_delaunay<3>(numPoints, numOutputSteps, outputVTKFile);
  }
  else
  {
    params.showhelp();
  }

  return EXIT_SUCCESS;
}
