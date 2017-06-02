/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef MINT_FEM_BASIS_TYPES_HPP_
#define MINT_FEM_BASIS_TYPES_HPP_

/*!
 *******************************************************************************
 * \enum FEBasisTypes
 *
 * \brief Enumerates the different types of supported Finite Element basis.
 *
 * \see FEBasis
 *******************************************************************************
 */
enum FEBasisTypes {
  MINT_UNDEFINED_BASIS=-1, /*!< Undefined basis type */
  MINT_LAGRANGE_BASIS,     /*!< Lagrange basis type */

  MINT_NUM_BASIS           /*!< MINT_NUM_BASIS */
};

namespace axom {
namespace mint {

/*!
 *******************************************************************************
 * \brief Array of strings corresponding to each Finite Element Basis.
 *
 * \note The length of the array is MINT_NUM_BASIS
 *
 * \note It is used to get a string representation for a Finite Element basis,
 *  e.g., for, debugging etc.
 *******************************************************************************
 */
static const char* basis_name[] = {
   "LAGRANGE_BASIS",
};


}
}



#endif
