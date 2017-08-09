// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file delaunay_triangulation.cpp
 * \brief This file contains an example that creates Delaunay triangulation
 * sequentially from a set of random points contained within a rectangle
 * centered at the origin.
 * It uses SLAM IA, a topological mesh data structure.
 *******************************************************************************
 */

#include "axom/slic.hpp"
#include "axom/quest/Delaunay.hpp"

#include <time.h>
#include <iomanip>

using namespace axom;

double BOUNDING_BOX_DIM[] = {
  2.0,
  2.0,
  2.0};  //This is the length of the bounding rectangle to contain the random points

typedef int IndexType;
typedef double DataType;

typedef axom::primal::Point2D Point2D;
typedef axom::primal::Point3D Point3D;

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
    std::cout << "Argument usage:" << std::endl
              << "  --help(-h)             Show this help message." << std::endl
              << "  --nrandpt(-n) [N]      The number of random points.  "
                 "Default N = 10.  "
              << std::endl
              << "  --dim(-d) [D]          The dimension of the mesh. 2 = "
                 "triangles in 2D, 3 = tetrahedrons in 3D."
              << std::endl
              << "  --outputsteps(-s)      The intermediate steps will be "
                 "written out as VTK files"
              << std::endl
              << "  --outfile(-o) [fname]  The VTK output file. Default fname "
                 "= delaunay_output."
              << std::endl
              << std::endl;
  };
};

Input::Input(int argc, char** argv)
  : outputVTKFile("delaunay_output")
  , numRandPoints(10)
  , outputSteps(false)
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

  SLIC_INFO("Using parameter values: "
            << std::endl
            << "  nrandpt = " << numRandPoints << std::endl
            << "  dimension = " << dimension << std::endl
            << "  outputSteps = " << outputSteps << std::endl
            << "  outfile = " << outputVTKFile << std::endl);
}

template <unsigned int dimension>
int run_delaunay(int num_points, bool outputVTKsteps, std::string outputVTKFile)
{
  // Create Delaunay Triangulation on random points
  typedef axom::quest::Delaunay<dimension> Delaunay;
  typedef typename Delaunay::PointType PointType;
  typedef axom::primal::BoundingBox<DataType, dimension> BoundingBox;
  Delaunay dt;

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
  time_t rseed_val = time(NULL);
  SLIC_INFO("Using " << rseed_val << " as seed for rand()");
  srand(rseed_val);

  if(outputVTKsteps)
  {
    std::ostringstream stm;
    stm << outputVTKFile << "_box.vtk";
    dt.writeToVTKFile(stm.str());
  }

  for(int pt_i = 0; pt_i < num_points; pt_i++)
  {
    //SLIC_INFO("Adding point #" << pt_i );
    //double x = (double)rand() / (double)RAND_MAX * BOUNDING_BOX_DIM[0] - BOUNDING_BOX_DIM[0]/2.0;
    //double y = (double)rand() / (double)RAND_MAX * BOUNDING_BOX_DIM[1] - BOUNDING_BOX_DIM[1]/2.0;

    double pt_coord[dimension];
    for(unsigned int d = 0; d < dimension; d++)
      pt_coord[d] = (double)rand() / (double)RAND_MAX * BOUNDING_BOX_DIM[d] -
        BOUNDING_BOX_DIM[d] / 2.0;
    PointType new_pt(pt_coord, dimension);

    dt.insertPoint(new_pt);

    if(outputVTKsteps)
    {
      std::ostringstream stm;
      stm << outputVTKFile << std::setfill('0') << std::setw(2) << pt_i << ".vtk";
      dt.writeToVTKFile(stm.str());
    }
  }

  dt.removeBoundary();  //Remove the starting rectangular box

  // Write the final mesh to a vtk file
  std::ostringstream stm;
  stm << outputVTKFile << ".vtk";
  dt.writeToVTKFile(stm.str());

  //dt.printMesh();

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

  const bool outputVTKsteps = params.outputSteps;
  const int num_points = params.numRandPoints;
  const std::string outputVTKFile = params.outputVTKFile;
  const int dimension = params.dimension;

  if(dimension == 2)
  {
    run_delaunay<2>(num_points, outputVTKsteps, outputVTKFile);
  }
  else if(dimension == 3)
  {
    run_delaunay<3>(num_points, outputVTKsteps, outputVTKFile);
  }
  else
  {
    params.showhelp();
  }

  return EXIT_SUCCESS;
}
