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

#ifndef MINT_FEM_BASIS_TYPES_HPP_
#define MINT_FEM_BASIS_TYPES_HPP_

#include <string>

/*!
 * \enum FEBasisTypes
 *
 * \brief Enumerates the different types of supported Finite Element basis.
 *
 * \see FEBasis
 */
enum FEBasisTypes
{
  MINT_UNDEFINED_BASIS=-1, /*!< Undefined basis type */
  MINT_LAGRANGE_BASIS,     /*!< Lagrange basis type */

  MINT_NUM_BASIS_TYPES           /*!< MINT_NUM_BASIS */
};

namespace axom
{
namespace mint
{

/*!
 * \brief Array of strings corresponding to each Finite Element Basis.
 *
 * \note The length of the array is MINT_NUM_BASIS
 *
 * \note It is used to get a string representation for a Finite Element basis,
 *  e.g., for, debugging etc.
 */
static const std::string basis_name[] = {
  "LAGRANGE_BASIS",
};

}
}

#endif
