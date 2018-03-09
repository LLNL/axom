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

#ifndef MINT_SRC_UTILS_VTK_UTILS_HPP
#define MINT_SRC_UTILS_VTK_UTILS_HPP

#include <string> // for std::string

namespace axom
{
namespace mint
{

class Mesh;

/*!
 * \brief Writes a mesh to a VTK file using the legacy ASCII format.
 *  that can be visualized with VisIt or ParaView.
 * \param [in] mesh the mesh to write out.
 * \param [in] file_path the path of the file to write to.
 * \return an error code, zero signifies a successful write.
 * \pre mesh != AXOM_NULLPTR
 * \note Thise method is primarily intended for debugging.
 */
int write_vtk( const Mesh* mesh, const std::string& file_path );

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_SRC_UTILS_VTK_UTILS_HPP */
