/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file Triangle.hpp
 *
 * \date Aug 19, 2015
 *******************************************************************************
 */


// ATK Toolkit includes
// #include "common/ATKMacros.hpp"
#include "axom/Types.hpp"
#include "axom_utils/FileUtilities.hpp"
#include "axom_utils/Timer.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/BVHTree.hpp"
#include "primal/HyperSphere.hpp"
#include "primal/orientation.hpp"
#include "primal/Point.hpp"
#include "quest/STLReader.hpp"
// #include "quest/SquaredDistance.hpp"  // <<!!???   
#include "primal/Triangle.hpp"
#include "primal/Vector.hpp"
#include "quest/SignedDistance.hpp"
#include "primal/intersection.hpp"
#include "primal/UniformGrid.hpp"

#include "mint/UnstructuredMesh.hpp"
#include "mint/UniformMesh.hpp"
#include "mint/Mesh.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/FieldData.hpp"
#include "mint/Field.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"

#include "slam/Utilities.hpp"


// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdio.h>

using namespace axom;

typedef mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
typedef primal::Triangle<double, 3> Triangle3;

typedef primal::Point<double, 3> Point3;
typedef primal::BoundingBox<double, 3> SpatialBoundingBox;
typedef primal::UniformGrid<int, 3> UniformGrid3;
typedef primal::Vector<double, 3> Vector3;
typedef primal::Segment<double, 3> Segment3;


typedef struct Input{
  std::string stlInput;
  std::string textOutputFile;
  int resolution;
}Input;

SpatialBoundingBox compute_bounds( mint::Mesh* mesh);
Triangle3 getMeshTriangle(int i,  mint::Mesh* surface_mesh);
bool shareVertices(Triangle3 t1, Triangle3 t2);
bool trianglesAreCoplanar(Triangle3 t1, Triangle3 t2);
bool trianglesOverlap(Triangle3 t1, Triangle3 t2);
bool checkTT(Triangle3& t1, Triangle3& t2, std::ofstream& ofs, int i, int j);
void naiveAlgorithm(std::ofstream& ofs, mint::Mesh* surface_mesh);
void vgAlgorithm(std::ofstream& ofs, mint::Mesh* surface_mesh, int resolution);
void init_params(Input& params,int argc, char ** argv);
bool trianglesIn(Triangle3 t1, Triangle3 t2);
inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS);

//compute_bounds is from code reused from George Zagaris' sphere example.  Computes the bounding box of a mesh
SpatialBoundingBox compute_bounds( mint::Mesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  SpatialBoundingBox meshBB;
  Point3 pt;

  for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
  {
    mesh->getMeshNode( i, pt.data() );
    meshBB.addPoint( pt );
  } // END for all nodes

  SLIC_ASSERT( meshBB.isValid() );

  return meshBB;
}

//helper function that retrieves the triangle at the given index
Triangle3 getMeshTriangle(int i,  mint::Mesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getMeshNumberOfCellNodes(i) == 3);
  primal::Point<int, 3> triCell;
  surface_mesh->getMeshCell( i, triCell.data());
  primal::Point< double,3 > A1;
  primal::Point< double,3 > B1;
  primal::Point< double,3 > C1;

  surface_mesh->getMeshNode( triCell[0],A1.data() );
  //SLIC_INFO("A is " <<A1);
  surface_mesh->getMeshNode( triCell[1],B1.data() );
  //SLIC_INFO("B is " <<B1);
  surface_mesh->getMeshNode( triCell[2],C1.data() );
  //SLIC_INFO("C is " <<C1);
  Triangle3 t1 = Triangle3(A1,B1,C1);
  //SLIC_INFO("T is " <<t1);

  return t1;
}
inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS=1.0e-9) {
  return axom::utilities::isNearlyEqual(p1[0], p2[0], EPS) && \
    axom::utilities::isNearlyEqual(p1[1], p2[1], EPS) && \
    axom::utilities::isNearlyEqual(p1[2], p2[2], EPS);

}

//check if any of the vertices are the same point (all of their points are within 1.0e-11 of eachother)
bool shareVertices(Triangle3 t1, Triangle3 t2)
{
  //SLIC_INFO("t1 is "<< t1 << " t2 is "<< t2);
  return ((pointIsNearlyEqual(t1[0], t2[0])) || \
          (pointIsNearlyEqual(t1[0], t2[1])) || \
          (pointIsNearlyEqual(t1[0], t2[2])) || \
          (pointIsNearlyEqual(t1[1], t2[0])) || \
          (pointIsNearlyEqual(t1[1], t2[1])) || \
          (pointIsNearlyEqual(t1[1], t2[2])) || \
          (pointIsNearlyEqual(t1[2], t2[0])) || \
          (pointIsNearlyEqual(t1[2], t2[1])) || \
          (pointIsNearlyEqual(t1[2], t2[2])));

}

//It turns out that I didn't need the below function... I'm leaving this here in case it is needed after all
bool trianglesAreCoplanar(Triangle3 t1, Triangle3 t2) {
  //Step 11.1: Check if the triangles are coplanar
  Vector3 t1Normal = Vector3::cross_product(Vector3(t1[0], t1[1]), Vector3(t1[0], t1[2]));
  Vector3 t2Normal = Vector3::cross_product(Vector3(t2[2], t2[0]), Vector3(t2[2], t2[1]));

  double dAi = (Vector3(t2[2], t1[0])).dot(t2Normal);
  double dBi = (Vector3(t2[2],t1[1])).dot(t2Normal);
  double dCi = (Vector3(t2[2],t1[2])).dot(t2Normal); 

  double dAj = (Vector3(t1[2],t2[0])).dot(t1Normal);
  double dBj = (Vector3(t1[2],t2[1])).dot(t1Normal);
  double dCj = (Vector3(t1[2],t2[2])).dot(t1Normal);
    
  return (axom::utilities::isNearlyEqual(dAi, 0.0, 1.0e-12) && \
          axom::utilities::isNearlyEqual(dBi, 0.0, 1.0e-12) && \
          axom::utilities::isNearlyEqual(dCi, 0.0, 1.0e-12)) || \
    (axom::utilities::isNearlyEqual(dAj, 0.0, 1.0e-12) && \
     axom::utilities::isNearlyEqual(dBj, 0.0, 1.0e-12) && \
     axom::utilities::isNearlyEqual(dCj, 0.0, 1.0e-12)); 
}

//Requires t1 and t2 share at least 1 vertex
bool trianglesIn(Triangle3 t1, Triangle3 t2) {
  /*Step 11.1: given that t1 and t2 share a vertex, permute all possible vertex matching pairs
    We have three cases: 
    (1) triangles share all 3 vertices (are the same triangle)
    (2) triangles share 2 vertices-- if they overlap not on the edge formed by these two vertices, 
    it is because one of their vertices lies within the other. Note that this function must be called twice to check this case, as only one of 
    the non-intersecting vertices need lie in the other triangle for an intersection to occur. 
    (3) Triangles share 1 vertex -- if they overlap or intersect, one of the triangles two non-intersecting vertices must form
    a segment that intersects the other. Note that this function must be called twice to check this case, as only one of the segments
    formed by the non-intersecting vertices need intersect the other triangle for an intersection to occur.
  */
  bool returning;
  if ((pointIsNearlyEqual(t1[0],t2[0])) || (pointIsNearlyEqual(t1[0],t2[1])) ||(pointIsNearlyEqual(t1[0],t2[2]))) {
    if ((pointIsNearlyEqual(t1[1],t2[0])) || (pointIsNearlyEqual(t1[1],t2[1])) ||(pointIsNearlyEqual(t1[1],t2[2]))){
      if ((pointIsNearlyEqual(t1[2],t2[0])) || (pointIsNearlyEqual(t1[2],t2[1])) ||(pointIsNearlyEqual(t1[2],t2[2]))) {
        return true; //Case 1: triangles are the same triangle, so they must be intersecting in a region
      }
      else {
        returning= trianglesAreCoplanar(t1, t2) && (t2.checkInTriangle(t1[2]));  //Case 2: t1.B and t1.A are a side on t2
        return returning;
      }
        
    }
    else if ((pointIsNearlyEqual(t1[2],t2[0])) || (pointIsNearlyEqual(t1[2],t2[1])) ||(pointIsNearlyEqual(t1[2],t2[2]))) {
      returning=  trianglesAreCoplanar(t1, t2) && (t2.checkInTriangle(t1[1])); //t1.C and t1.A are a side on t2
      return returning;
    }
    else {
      Segment3 s= Segment3(t1[1], t1[2]);
      returning= primal::intersect(t2, s);  //t1.A is a shared vertex
      return returning;
    }

  }
  else if ((pointIsNearlyEqual(t1[1],t2[0])) || (pointIsNearlyEqual(t1[1],t2[1])) ||(pointIsNearlyEqual(t1[1],t2[2]))) {  //t1.A is not a shared vertex, at least one of t1.B or t1.C share a vertex with t2
    if ((pointIsNearlyEqual(t1[2],t2[0])) || (pointIsNearlyEqual(t1[2],t2[1])) ||(pointIsNearlyEqual(t1[2],t2[2]))){
      returning= trianglesAreCoplanar(t1, t2) && (t2.checkInTriangle(t1[0]));  //t1.B and t1.C are a side on t2
      return returning;
    }
    else {
      Segment3 s= Segment3(t1[0], t1[2]);
      returning= primal::intersect(t2, s);  //t1.B is a shared vertex,  
      return returning;
    }
  }
  else  { //t1.C must be a shared vertex, because t1.A and t1.B are not shared
    Segment3 s= Segment3(t1[0], t1[1]);
    returning= primal::intersect(t2, s);  //t1.B is a shared vertex,  
    return returning;
  }
}

//Requires t1 and t2 share at least 1 vertex
bool trianglesOverlap(Triangle3 t1, Triangle3 t2) {
  return ((trianglesIn(t1, t2) || trianglesIn(t2,t1)));
}

bool checkTT(Triangle3& t1, Triangle3& t2, std::ofstream& ofs, int i, int j)
{
  //Step 9:  Ignore these degenerate triangles, as we will eventually write out all
  //degenerate triangles in main for loop
  if (t2.degenerate()) return false;
  //Step 10:  If two triangles share at least one vertice, they are either connected only along an edge/vertex or one of the 
  //triangles has at least one vertex lying on the interior or along an edge (excluding the vertices) of the other triangle
  if (shareVertices(t1, t2)) {
    if (trianglesOverlap(t1, t2)) {
      ofs<<"\n INTERSECTION FOUND Triangle 1 at index "<< i<<" is "<<t1<<"and triangle 2 at index "<<j<<"is "<< t2 << "\n";
      return true;
    }
    return false;
  } 
  //Step 11.2: If no shared vertice, run the standard intersection test
  return  primal::intersect(t1,  t2);
}

void naiveAlgorithm(std::ofstream& ofs, mint::Mesh* surface_mesh) {
        


  //Step 6-11 overview: For each triangle, check for intersections against every other 
  //triangle with a greater index in the mesh, excluding degenerate triangles and
  //triangles whose only intersection with the initial triangle occurs at a vertice
  const int ncells = surface_mesh->getMeshNumberOfCells();
  SLIC_INFO("Checking mesh with a total of "<< ncells<< " cells.");

  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();
  int counter=0; //total intersections found, will be written to outfile at end
  int degenerateCounter=0;  //total degenerate triangles found
  //Step 6:  Iterate over all the triangles in the mesh
  for (int i=0; i< ncells; i++) {
    t1=getMeshTriangle(i,  surface_mesh);
    //Step 7: Check if the triangle is degenerate, if so, write out to output file and skip
    if (t1.degenerate()) {
      degenerateCounter++;     
      ofs<<"Warning, degenerate triangle at index "<< i<<" with value "<<t1;
      continue;
    }
    //Step 8:  If the triangle is not degenerate, test against all other triangles that this
    //triangle has not been checked against
    for (int j=i+1; j< ncells; j++) {
      t2=getMeshTriangle(j,  surface_mesh);
      //Step into CheckTT for steps 9-11
      if (checkTT(t1, t2, ofs, i, j)) counter++;
    }
  }
  //step 12, write results to outfile
  ofs<<"A total of "<< counter <<  "  triangle intersections were found.";
  ofs<<"A total of "<< degenerateCounter <<  " degenerate triangles were found.";
}

void vgAlgorithm(std::ofstream& ofs, mint::Mesh* surface_mesh, int resolution) {


  int counter=0;  //total intersections found, will be written to outfile at end
  int degenerateCounter=0; //total degenerate triangles found
  //Step 6:  Construct our virtual grid
  int intersectingTriangleCount=0;
  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();
  SLIC_INFO("Running mesh_tester with virtual_grid");

  //Create a bounding box around mesh to find the minimum point
  SpatialBoundingBox meshBB  = compute_bounds(surface_mesh);
  const Point3 minBBPt= meshBB.getMin();
  const Point3 maxBBPt= meshBB.getMax();

  //figure out the appropriate side sizes for our given resolution
  const double sideSizes [3] = {((maxBBPt[0]-minBBPt[0]+1.0e-12)/((double)resolution)), \
                                ((maxBBPt[1]-minBBPt[1]+1.0e-12)/((double)resolution)), \
                                ((maxBBPt[2]-minBBPt[2]+1.0e-12)/((double)resolution))};

  //create an array containing the resolution for each side.  This will be uniform for now
  //WIP client ability to specify x,y and z resolutions
  int resolutions[3]={resolution,resolution,resolution};

  std::cerr<<"Building virtual grid...\n";
  UniformGrid3 vg(minBBPt.data(), maxBBPt.data(), resolutions);
  const int ncells = surface_mesh->getMeshNumberOfCells();
  for (int i=0; i< ncells; i++) {
    if (i%(ncells/100)==0) {
      std::cerr<<"Building grid is "<<100.0*(double(i)/double(ncells)) <<" percent done \n";

    }
    SpatialBoundingBox triBB;
    t1=getMeshTriangle(i,  surface_mesh);
    triBB.addPoint(t1[0]);
    triBB.addPoint(t1[1]);
    triBB.addPoint(t1[2]);

    vg.insert(triBB, i);
  }


  //Step 7:  Iterate through triangle indices from first index to last index and check against any other triangles with indexes 
  //greater than the index z who also share a virtual grid bin with the triangle at index z.  Avoid repeat checks by
  //only checking against triangles with indexes greater than z. 
  SLIC_INFO("Checking mesh with a total of "<< ncells<< " cells.");
  for (size_t z=0; z< ncells; z++) {
    if (z%(ncells/100)==0) {
      //the below will not actually reflect the time the algorithm will take perfectly, as triangles at earlier indexes
      //must do more checks than triangles at later indexes
      std::cerr<<"Querying grid is "<<100.0*(double(z)/double(ncells)) <<" percent done \n";
    }
        
    //Step 7.1:  Retrieve the triangle at index z and construct a bounding box around it
    SpatialBoundingBox triBB2;
    t1=getMeshTriangle(z,  surface_mesh);
    //Step 7.2:  Check if this triangle is degenerate, if so, write to ofs and skip this triangle test
    if (t1.degenerate()) { 
      degenerateCounter++;    
      ofs<<"Warning, degenerate triangle at index "<< z<<" with ~value "<<t1;
      continue;
    }
    triBB2.addPoint(t1[0]);
    triBB2.addPoint(t1[1]);
    triBB2.addPoint(t1[2]);

    Point3 minBBPt2,maxBBPt2;
        
    minBBPt2 = triBB2.getMin();
    maxBBPt2 = triBB2.getMax();   
    bool notBadBefore=true;
    //step 7.3:  Find the virtual grid indices that this triangle will touch.  The below breaks the black box created by the virtual grid
    //This can be avoided by providing increased functinality in virtual grid
    std::vector<int> trianglesInBin;
    size_t lowerIndex = vg.getBinIndex(minBBPt2);
    size_t stepsPerSide [3] = {((size_t)(1.0+(maxBBPt2[0]-minBBPt2[0]+1.0e-12)/((double)sideSizes[0]))), \
                               ((size_t)(1.0+(maxBBPt2[1]-minBBPt2[1]+1.0e-12)/((double)sideSizes[1]))), \
                               ((size_t)(1.0+(maxBBPt2[2]-minBBPt2[2]+1.0e-12)/((double)sideSizes[2])))};
    //   SLIC_INFO("Steps per side is: "<< stepsPerSide[0] <<","<<stepsPerSide[1] <<","<<stepsPerSide[2] );

    //Step 7.4: iterate through all the bins that the triangle bounding box touches
    for (size_t i =0; i<stepsPerSide[0]; i++)  {
      for (size_t j =0; j<stepsPerSide[1]; j++) {
        for (size_t k =0; k<stepsPerSide[2]; k++) {
          //Step 7.5: retrieve triangles for the bin
          trianglesInBin= vg.getBinContents(lowerIndex + i + resolutions[0]*(j + (k * resolutions[1])));    
          int binSize=trianglesInBin.size();

          //Step 7.6: Iterate through triangles in the bin. Check if we have already tested this intection pair
          //by comparing the indices.  If we have, skip, otherwise, proceed to step 8
          for (size_t l=0; l<binSize; l++) {
            int t2Index=trianglesInBin[l];
            SLIC_ASSERT(trianglesInBin[l]>=0);
            if (t2Index<=z) continue;  
            else{ 
              //Step 8: if we haven't already tested the intersection, run the intersection test
              t2=getMeshTriangle(t2Index,  surface_mesh);
              //Step into CheckTT for steps 9-11
              if (checkTT(t1, t2, ofs, i, t2Index)) {
                counter++;
                if (notBadBefore) {
                  intersectingTriangleCount ++;
                  notBadBefore=false;
                }
              }
            }
          }
        }
      }
    }
  }
  //step 12, write results to outfile
  ofs<<"A total of "<< counter <<  " triangle intersections were found.";
  ofs<<"A total of "<< intersectingTriangleCount <<  " intersecting triangles were found. \n";
  ofs<<"A total of "<< degenerateCounter <<  " degenerate triangles were found.";
}

//helper function to update the parameters with the command line arguments
void init_params(Input& params,int argc, char ** argv){
  if((argc > 1)){
    std::string help = argv[1];
    if(help == "--help"){
      std::cout << "Possible flags/usage are as follows" << std::endl
                << "--resolution: resolution of virtual grid-- make this 1 if you want to run the naive algorithm" <<std::endl
                << "--infile -- the stl input file" << std::endl
                << "--outfile -- the text output file containing the results" << std::endl;
      exit(100);

    }
    if((argc % 2 ) != 1){
      std::cerr << "Incorrect number of command line arguments, needs to be odd";
      exit(200);
    }
    for(int i = 1; i< argc;){
      std::string arg = argv[i];
      if(arg == "--resolution"){
        params.resolution = atoi(argv[++i]);
      }
      else if(arg == "--infile"){
        params.stlInput = argv[++i];
      }
      else if(arg == "--outfile"){
        params.textOutputFile = argv[++i];
      }
      ++i;
    }
  }
  std::cerr << "Parameters are: " << std::endl
            << "resolution = " << params.resolution << std::endl
            << "infile = " << params.stlInput << std::endl
            << "outfile = " << params.textOutputFile << std::endl;
}




int main( int argc, char** argv ) {
  // STEP 0: Initialize SLIC Environment
  slic::initialize();
  slic::setLoggingMsgLevel( axom::slic::message::Debug );

  // Create a more verbose message for this application (only level and message)
  std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";
  slic::GenericOutputStream* defaultStream =
    new slic::GenericOutputStream(&std::cout);
  slic::GenericOutputStream* compactStream =
    new slic::GenericOutputStream(&std::cout, slicFormatStr);
  slic::addStreamToMsgLevel(defaultStream, axom::slic::message::Error);
  slic::addStreamToMsgLevel(compactStream, axom::slic::message::Warning);
  slic::addStreamToMsgLevel(compactStream, axom::slic::message::Info);
  slic::addStreamToMsgLevel(compactStream, axom::slic::message::Debug);

  //Step 1: Initialize default parameters and update parameters with command line arguments:
  Input params;
  const std::string defaultFileName = "plane.stl";
  const std::string defaultDir = "/g/g17/edesanto/axom/src/components/primal/data/";
  params.stlInput = utilities::filesystem::joinPath(defaultDir, defaultFileName);
  params.textOutputFile="meshTestResults.txt";
  params.resolution=10;
  init_params(params,argc,argv);
  SLIC_ASSERT(params.resolution>0);
   

  // STEP 2: read file
  SLIC_INFO("Reading file: " <<  params.stlInput << "...\n");
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName(  params.stlInput );
  reader->read();
  SLIC_INFO("done\n");

  // STEP 3: get surface mesh
  mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader-> getMesh( static_cast<TriangleMesh*>( surface_mesh ) );

  // STEP 4: Delete the reader
  delete reader;
  reader = AXOM_NULLPTR;

  // STEP 5: write vtk file
  //write_vtk( surface_mesh, "surface_mesh.vtk" );
  std::ofstream ofs;
  ofs.open(params.textOutputFile);

  if (params.resolution==1) {
    // Naive method -- check inside for steps 5-12
    naiveAlgorithm(ofs, surface_mesh);
  }
  else {
    //  Virtual grid method -- check inside for steps 5-12
    vgAlgorithm(ofs, surface_mesh, params.resolution);
  }

  ofs.close();

}
