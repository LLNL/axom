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

#ifndef MINT_MESHTYPES_HPP_
#define MINT_MESHTYPES_HPP_

namespace axom
{
namespace mint
{

/*!
 * \brief Defines the basic mesh types supported by mint.
 */
enum MeshTypes
{
  UNDEFINED_MESH = -1,   ///< UNDEFINED_MESH

  UNSTRUCTURED_MESH,     ///< UNSTRUCTURED_MESH

  STRUCTURED_MESH,       ///< STRUCTURED_MESH
  RECTILINEAR_MESH,      ///< RECTILINEAR_MESH
  UNIFORM_MESH,          ///< UNIFORM_MESH
  PARTICLE_MESH,         ///< PARTICLE_MESH

  NUM_MESH_TYPES         ///< NUM_MESH_TYPES
};

} /* namespace mint */
} /* namespace axom */
#endif /* MINT_MESHTYPE_HPP_ */
