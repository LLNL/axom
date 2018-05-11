/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*! \file primal_introduction.cpp
 *  \brief This example code is a basic demonstration of how to use
 *  Primal to represent geometric primitives, perform geometric operations,
 *  and use a spatial index.
 */

/* This example code contains snippets used in the Primal Sphinx documentation.
 * They begin and end with comments
 *
 * first_example_header_start
 * first_example_header_end
 * first_example_using_start
 * first_example_using_end
 * first_example_primitive_start
 * first_example_primitive_end
 * first_example_operation_start
 * first_example_operation_end
 * first_example_naive_intersection_start
 * first_example_naive_intersection_end
 * first_example_ugrid_intersection_start
 * first_example_ugrid_intersection_end
 *
 * each prepended with an underscore.
 */

// _first_example_header_start
// Axom primitives
#include "primal/BoudingBox.hpp"
#include "primal/Plane.hpp"
#include "primal/Point.hpp"
#include "primal/Polygon.hpp"
#include "primal/Ray.hpp"
#include "primal/Triangle.hpp"

// Axom operations
#include "primal/clip.hpp"
#include "primal/closest_point.hpp"
#include "primal/compute_bounding_box.hpp"
#include "primal/intersect.hpp"
#include "primal/orientation.hpp"
#include "primal/squared_distance.hpp"

// Axom spatial index
#include "primal/UniformGrid.hpp"
// _first_example_header_end

// C++ headers
#include <cmath> // do we need this?

// _first_example_using_start
// "using" directives to simplify code
using namespace axom;
using namespace primal;
// _first_example_using_end

DataStore* create_datastore(int* region)
{
}

int main(int argc, char** argv)
{

  // Deal with unused variables
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);

  int region[3375];

  DataStore* ds = create_datastore(region);
  access_datastore(ds);

  DataStore* tds = create_tiny_datastore();
  save_as_blueprint(tds);

  return 0;
}
