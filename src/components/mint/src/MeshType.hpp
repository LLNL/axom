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
#ifndef MESHTYPE_HPP_
#define MESHTYPE_HPP_


namespace mint {

/*!
 * \brief Defines the basic mesh types supported by mint
 */
enum MeshType {
  UNDEFINED_MESH = -1,             //!< UNDEFINED_MESH

  UNSTRUCTURED_SEGMENT_MESH,       //!< UNSTRUCTURED_SEGMENT_MESH
  UNSTRUCTURED_ALL_TRIANGLE_MESH,  //!< UNSTRUCTURED_TRIANGLE_MESH
  UNSTRUCTURED_ALL_QUAD_MESH,      //!< UNSTRUCTURED_QUAD_MESH
  UNSTRUCTURED_ALL_TET_MESH,       //!< UNSTRUCTURED_TET_MESH
  UNSTRUCTURED_ALL_HEX_MESH,       //!< UNSTRUCTURED_HEX_MESH
  UNSTRUCTURED_MIXED_ELEMENT_MESH, //!< UNSTRUCTURED_MIXED_ELEMENT_MESH

  STRUCTURED_CURVILINEAR_MESH,     //!< STRUCTURED_CURVILINEAR_MESH
  STRUCTURED_RECTILINEAR_MESH,     //!< STRUCTURED_RECTILINEAR_MESH
  STRUCTURED_UNIFORM_MESH,         //!< STRUCTURED_UNIFORM_MESH

  PARTICLE_MESH,                   //!< PARTICLE_MESH

  NUM_MESH_TYPES                   //!< NUM_MESH_TYPES
};

namespace mesh_properties {

/*!
 *******************************************************************************
 * \brief Defines corresponding mesh type for each cell type.
 * \note The array is mint::NUM_CELL_TYPES long
 * \see mint::CellType
 *******************************************************************************
 */
static const int mesh_of_cell_type[] = {

    PARTICLE_MESH,                   // VERTEX,

    UNSTRUCTURED_SEGMENT_MESH,       // LINE_SEGMENT
    UNSTRUCTURED_ALL_TRIANGLE_MESH,  // LINEAR_TRIANGLE,
    UNSTRUCTURED_ALL_QUAD_MESH,      // LINEAR_QUAD,
    UNSTRUCTURED_ALL_TET_MESH,       // LINEAR_TET,
    UNSTRUCTURED_ALL_HEX_MESH,       // LINEAR_HEX,

    UNDEFINED_MESH,                  // LINEAR_PRISM,
    UNDEFINED_MESH,                  // LINEAR_PYRAMID,

    UNSTRUCTURED_MIXED_ELEMENT_MESH, // MIXED
};

/*!
 *******************************************************************************
 * \brief Defines whether a mesh type stores explicit coordinates.
 * \note The array is mint::NUM_MESH_TYPES long
 * \see mint::MeshType
 *******************************************************************************
 */
static const bool explicit_coordinates[] = {
    true,  // UNSTRUCTURED_SEGMENT_MESH
    true,  // UNSTRUCTURED_ALL_TRIANGLE_MESH
    true,  // UNSTRUCTURED_ALL_QUAD_MESH
    true,  // UNSTRUCTURED_ALL_TET_MESH
    true,  // UNSTRUCTURED_ALL_HEX_MESH
    true,  // UNSTRUCTURED_MIXED_ELEMENT_MESH

    true,  // STRUCTURED_MESH
    true,  // RECTILINEAR_MESH
    false, // UNIFORM_MESH

    true,  // PARTICLE_MESH
};

/*!
 *******************************************************************************
 * \brief Defines whether a mesh type stores explicit connectivity.
 * \note The array is mint::NUM_MESH_TYPES long
 * \see mint::MeshType
 *******************************************************************************
 */
static const bool explicit_connectivity[] = {
    true,  // UNSTRUCTURED_SEGMENT_MESH
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


} /* namespace mesh_properties */

} /* namespace mint */


#endif /* MESHTYPE_HXX_ */
