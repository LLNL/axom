/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

#ifndef MESHTYPE_HPP_
#define MESHTYPE_HPP_

/*!
 * \brief Defines the basic mesh types supported by mint.
 */
enum
{
  MINT_UNDEFINED_MESH = -1,             //!< UNDEFINED_MESH

  MINT_UNSTRUCTURED_SEGMENT_MESH,       //!< UNSTRUCTURED_SEGMENT_MESH
  MINT_UNSTRUCTURED_TRIANGLE_MESH,      //!< UNSTRUCTURED_TRIANGLE_MESH
  MINT_UNSTRUCTURED_QUAD_MESH,          //!< UNSTRUCTURED_QUAD_MESH
  MINT_UNSTRUCTURED_TET_MESH,           //!< UNSTRUCTURED_TET_MESH
  MINT_UNSTRUCTURED_HEX_MESH,           //!< UNSTRUCTURED_HEX_MESH
  MINT_UNSTRUCTURED_MIXED_ELEMENT_MESH, //!< UNSTRUCTURED_MIXED_ELEMENT_MESH

  MINT_STRUCTURED_CURVILINEAR_MESH,     //!< STRUCTURED_CURVILINEAR_MESH
  MINT_STRUCTURED_RECTILINEAR_MESH,     //!< STRUCTURED_RECTILINEAR_MESH
  MINT_STRUCTURED_UNIFORM_MESH,         //!< STRUCTURED_UNIFORM_MESH

  MINT_PARTICLE_MESH,                   //!< PARTICLE_MESH

  MINT_NUM_MESH_TYPES                   //!< NUM_MESH_TYPES
};

namespace axom
{
namespace mint
{

namespace mesh_properties
{

/*!
 * \brief Defines corresponding mesh type for each cell type.
 * \note The array is mint::NUM_CELL_TYPES long
 * \see mint::CellType
 */
static const int mesh_of_cell_type[] = {

  MINT_PARTICLE_MESH,                 // VERTEX,

  MINT_UNSTRUCTURED_SEGMENT_MESH,     // LINE_SEGMENT
  MINT_UNSTRUCTURED_TRIANGLE_MESH,    // LINEAR_TRIANGLE,
  MINT_UNSTRUCTURED_QUAD_MESH,        // LINEAR_QUAD,
  MINT_UNSTRUCTURED_TET_MESH,         // LINEAR_TET,
  MINT_UNSTRUCTURED_HEX_MESH,         // LINEAR_HEX,

  MINT_UNDEFINED_MESH,                // LINEAR_PRISM,
  MINT_UNDEFINED_MESH,                // LINEAR_PYRAMID,

  MINT_UNSTRUCTURED_QUAD_MESH,        // Q2 quad mesh
  MINT_UNSTRUCTURED_HEX_MESH,         // Q2 hex mesh

  MINT_UNSTRUCTURED_MIXED_ELEMENT_MESH,   // MIXED
};

/*!
 * \brief Defines whether a mesh type stores explicit coordinates.
 * \note The array is mint::NUM_MESH_TYPES long
 * \see mint::MeshType
 */
static const bool explicit_coordinates[] = {
  true,    // UNSTRUCTURED_SEGMENT_MESH
  true,    // UNSTRUCTURED_ALL_TRIANGLE_MESH
  true,    // UNSTRUCTURED_ALL_QUAD_MESH
  true,    // UNSTRUCTURED_ALL_TET_MESH
  true,    // UNSTRUCTURED_ALL_HEX_MESH
  true,    // UNSTRUCTURED_MIXED_ELEMENT_MESH

  true,    // STRUCTURED_MESH
  true,    // RECTILINEAR_MESH
  false,   // UNIFORM_MESH

  true,    // PARTICLE_MESH
};

/*!
 * \brief Defines whether a mesh type stores explicit connectivity.
 * \note The array is mint::NUM_MESH_TYPES long
 * \see mint::MeshType
 */
static const bool explicit_connectivity[] = {
  true,    // UNSTRUCTURED_SEGMENT_MESH
  true,    // UNSTRUCTURED_ALL_TRIANGLE_MESH
  true,    // UNSTRUCTURED_ALL_QUAD_MESH
  true,    // UNSTRUCTURED_ALL_TET_MESH
  true,    // UNSTRUCTURED_ALL_HEX_MESH
  true,    // UNSTRUCTURED_MIXED_ELEMENT_MESH

  false,    // STRUCTURED_MESH
  false,    // RECTILINEAR_MESH
  false,    // UNIFORM_MESH

  false,    // PARTICLE_MESH
};

} /* namespace mesh_properties */
} /* namespace mint */
} /* namespace axom */

#endif /* MESHTYPE_HPP_ */
