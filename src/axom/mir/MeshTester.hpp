// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MeshTester.hpp
 * 
 * \brief Contains the specification for the MeshTester class.
 * 
 */

#ifndef __MESH_TESTER_H__
#define __MESH_TESTER_H__

#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"  // unified header for slam classes and functions

#include "MIRMesh.hpp"

#include <algorithm>

namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{
  /**
   * \class MeshTester
   * 
   * \brief A class used to generate MIRMeshs with specific properties so
   *        that the reconstruction output can be validated visually.
   * 
   */
  class MeshTester
  {
    public:
      /**
      * \brief Default constructor. 
      */
      MeshTester();

      /**
       * \brief Default destructor.
       */
      ~MeshTester();

    public:
      /**
       * \brief Initializes an MIRMesh based on the example from Meredith 2004 paper.
       * 
       * \note The mesh is a 3x3 uniform grid of quads with 2 materials.
       * 
       * \return  The generated mesh.
       */
      MIRMesh  initTestCaseOne();

      /**
       * \brief Initializes an MIRMesh based on the example from Meredith and Childs 2010 paper.
       * 
       * \note The mesh is a 3x3 uniform grid of quads with 3 materials.
       * 
       * \return  The generated mesh.
       */
      mir::MIRMesh  initTestCaseTwo();

      /**
       * \brief Initializes an MIRMesh used for testing triangle clipping cases.
       * 
       * \note The mesh is a set of four triangles with 2 materials
       * 
       * \return  The generated mesh.
       */
      mir::MIRMesh  initTestCaseThree();

      /**
       * \brief Intializes a mesh used for testing a single circle of one materials surrounded by another.
       * 
       * \note The mesh is a 3x3 uniform grid with 2 materials and has a single circle in the center composed
       *       of one material and is surrounded by a second material.
       * 
       * \return  The generated mesh.
       */
      mir::MIRMesh  initTestCaseFour();

      /**
       * \brief Initializes a mesh to be used for testing a set of concentric circles centered in a uniform grid.
       * 
       * \param gridSize  The number of elements in the width and the height of the uniform grid.
       * \param numCircles  The number of concentric circles that are centered in the grid.
       * 
       * \note  Each circle is composed of a different material.
       * 
       * \return  The generated mesh.
       */
      mir::MIRMesh  initTestCaseFive(int gridSize, int numCircles); // multiple materials, multiple concentric circles

      /**
       * \brief Initializes a mesh composed of a uniform grid with a circle of material in it.
       * 
       * \param gridSize  The number of elements in the width and height of the uniform grid.
       * \param circleCenter  The center point of the circle.
       * \param circleRadius  The radius of the circle.
       * 
       * \return  The generated mesh.
       */
      mir::MIRMesh  createUniformGridTestCaseMesh(int gridSize, mir::Point2 circleCenter, axom::float64 circleRadius);
      
      /**
       * \brief Intializes a mesh to be used for validating the results of quad clipping.
       * 
       * \note The mesh is a 3x3 uniform grid with 2 materials and element volume fraction such
       *       that the mesh would be split horizontally through the middle.
       * 
       * \return  The generated mesh.
       */
      mir::MIRMesh  initQuadClippingTestMesh();

    private:

      /**
       * \brief Calculate the distance between the two given points.
       * 
       * \param p0  The first point.
       * \param p1  The second point.
       * 
       * \return The distance between the two points.
       */
      axom::float64  distance(mir::Point2 p0, mir::Point2 p1);

      /**
       * \brief Generates a 2D uniform grid of n x n elements.
       * 
       * \param gridSize  The number of elements in the width and height of the uniform grid.
       */
      mir::CellData  generateGrid(int gridSize);

      /**
       * \brief Calculates the percent overlap between the given circle and quad.
       * 
       * \param gridSize  The size of the uniform grid which will be sampled over to check for overlap.
       * \param circleCenter  The center point of the circle.
       * \param circleRadius  The radius of the circle.
       * \param quadP0  The upper left vertex of the quad.
       * \param quadP1  The lower left vertex of the quad.
       * \param quadP2  The lower right vertex of the quad.
       * \param quadP3  The upper right vertex of the quad.
       * 
       * /return The percent value overlap of the circle and the quad between [0, 1].
       */
      axom::float64  calculatePercentOverlapMonteCarlo(int gridSize, mir::Point2 circleCenter, axom::float64 circleRadius, mir::Point2 quadP0, mir::Point2 quadP1, mir::Point2 quadP2, mir::Point2 quadP3);
      
      /**
       * \brief Calculates the number of corners of the quad that are within the circle.
       * 
       * \param circleCenter The center point of the circle.
       * \param circleRadius The radius of the circle.
       * \param quadP0  The upper left vertex of the quad.
       * \param quadP1  The lower left vertex of the quad.
       * \param quadP2  The lower right vertex of the quad.
       * \param quadP3  The upper right vertex of the quad.
       * 
       * \return The number of corners of the quad that are within the circle.
       */
      int  circleQuadCornersOverlaps(mir::Point2 circleCenter, axom::float64 circleRadius, mir::Point2 quadP0, mir::Point2 quadP1, mir::Point2 quadP2, mir::Point2 quadP3);
  };
}
}

#endif
