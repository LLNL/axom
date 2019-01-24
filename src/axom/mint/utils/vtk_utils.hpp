/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

#ifndef MINT_SRC_UTILS_VTK_UTILS_HPP
#define MINT_SRC_UTILS_VTK_UTILS_HPP

#include <string> // for std::string

namespace axom
{
namespace mint
{

// Forward Declarations
class Mesh;
class FiniteElement;

/*!
 * \brief Writes a mesh to a VTK file using the legacy ASCII format.
 *  that can be visualized with VisIt or ParaView.
 * \param [in] mesh the mesh to write out.
 * \param [in] file_path the path of the file to write to.
 * \return an error code, zero signifies a successful write.
 * \pre mesh != nullptr
 * \note This method is primarily intended for debugging.
 */
int write_vtk( const Mesh* mesh, const std::string& file_path );

/*!
 * \brief Writes a FiniteElement to a VTK file in the legacy ASCII format.
 *
 * \param [in] fe reference to the finite element object.
 * \param [in] file_path path to the file to write
 *
 * \note This method is primarily intended for debugging.
 *
 * \return rc return code, zero signifies success
 */
int write_vtk( mint::FiniteElement& fe, const std::string& file_path );

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_SRC_UTILS_VTK_UTILS_HPP */
