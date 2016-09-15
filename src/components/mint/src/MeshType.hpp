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
 * \file
 *
 * \date Sep 19, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */
#ifndef MESHTYPE_HXX_
#define MESHTYPE_HXX_


namespace mint {

enum MeshType {
  UNDEFINED_MESH = -1,

  UNSTRUCTURED_ALL_TRIANGLE_MESH,
  UNSTRUCTURED_ALL_QUAD_MESH,
  UNSTRUCTURED_ALL_TET_MESH,
  UNSTRUCTURED_ALL_HEX_MESH,
  UNSTRUCTURED_MIXED_ELEMENT_MESH,

  STRUCTURED_CURVILINEAR_MESH,
  STRUCTURED_RECTILINEAR_MESH,
  STRUCTURED_UNIFORM_MESH,

  PARTICLE_MESH,

  NUM_MESH_TYPES
};

namespace mesh_properties {

static const int mesh_of_cell_type[] = {

    PARTICLE_MESH,                  // VERTEX,
    UNDEFINED_MESH,                 // LINE (?)

    UNSTRUCTURED_ALL_TRIANGLE_MESH, // LINEAR_TRIANGLE,
    UNSTRUCTURED_ALL_QUAD_MESH,     // LINEAR_QUAD,
    UNSTRUCTURED_ALL_TET_MESH,      // LINEAR_TET,
    UNSTRUCTURED_ALL_HEX_MESH,      // LINEAR_HEX,

    UNDEFINED_MESH,                 // LINEAR_PRISM,
    UNDEFINED_MESH,                 // LINEAR_PYRAMID,

    UNSTRUCTURED_MIXED_ELEMENT_MESH, // MIXED
};

static const bool explicit_coordinates[] = {
    true,  // UNIFORM_MESH
    true,  // UNSTRUCTURED_ALL_QUAD_MESH
    true,  // UNSTRUCTURED_ALL_TET_MESH
    true,  // UNSTRUCTURED_ALL_HEX_MESH
    true,  // UNSTRUCTURED_MIXED_ELEMENT_MESH

    true,  // STRUCTURED_MESH
    true,  // RECTILINEAR_MESH
    false, // UNIFORM_MESH

    true,  // PARTICLE_MESH
};

static const bool explicit_connectivity[] = {
    true,  // UNSTRUCTURED_ALL_TRIANGLE_MESH
    true,  // UNSTRUCTURED_ALL_QUAD_MESH
    true,  // UNSTRUCTURED_ALL_TET_MESH
    true,  // UNSTRUCTURED_ALL_HEX_MESH
    true,  // UNSTRUCTURED_MIXED_ELEMENT_MESH

    false,  // STRUCTURED_MESH
    false,  // RECTILINEAR_MESH
    false,  // UNIFORM_MESH

    false,  // PARTICLE_MESH
};

static const bool mixed_celltype[] = {
    false,  // UNSTRUCTURED_ALL_TRIANGLE_MESH
    false,  // UNSTRUCTURED_ALL_QUAD_MESH
    false,  // UNSTRUCTURED_ALL_TET_MESH
    false,  // UNSTRUCTURED_ALL_HEX_MESH
    true,   // UNSTRUCTURED_MIXED_ELEMENT_MESH

    false,  // STRUCTURED_MESH
    false,  // RECTILINEAR_MESH
    false,  // UNIFORM_MESH

    false,  // PARTICLE_MESH
};

} /* namespace mesh_properties */

} /* namespace mint */


#endif /* MESHTYPE_HXX_ */
