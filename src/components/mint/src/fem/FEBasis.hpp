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

#ifndef MINT_FEM_BASIS_HPP_
#define MINT_FEM_BASIS_HPP_

// Mint includes
#include "mint/CellType.hpp"
#include "mint/Lagrange.hpp"
#include "mint/ShapeFunction.hpp"
#include "mint/FEBasisTypes.hpp"

/*!
 * \def REGISTER_LAGRANGE_BASIS( C )
 *
 * \brief Macro used to register a new Lagrange Finite Element. This macro
 *  expands to a specialization of the FEBasis trait class that binds a
 *  Finite Element basis to a cell type.
 */
#define REGISTER_LAGRANGE_BASIS( C )                                          \
  template < >                                                                  \
  struct FEBasis< MINT_LAGRANGE_BASIS, C > {                                    \
    static const int BasisType = MINT_LAGRANGE_BASIS;                             \
    typedef mint::ShapeFunction< mint::Lagrange< C > > ShapeFunctionType; \
  }

namespace axom
{
namespace mint
{

/*!
 * \brief FEBasis is a traits class that binds a Finite Element basis type,
 *  e.g., MINT_LAGRANGE_BASIS, to a particular cell type, e.g., MINT_QUAD.
 *
 * \note This is an empty class that is intended to be specialized
 *
 * \see ShapeFunction
 * \see FEBasisTypes
 */
template < int BasisType, int CellType >
struct FEBasis { };

// Lagrange Basis
REGISTER_LAGRANGE_BASIS(  MINT_QUAD );
REGISTER_LAGRANGE_BASIS(  MINT_TRIANGLE );
REGISTER_LAGRANGE_BASIS(  MINT_TET );
REGISTER_LAGRANGE_BASIS(  MINT_HEX );
REGISTER_LAGRANGE_BASIS(  MINT_PRISM );
REGISTER_LAGRANGE_BASIS(  MINT_PYRAMID );

REGISTER_LAGRANGE_BASIS(  MINT_QUAD9 );
REGISTER_LAGRANGE_BASIS(  MINT_HEX27 );

} /* namespace mint */
} /* namespace axom */

// undef internal macros
#undef REGISTER_LAGRANGE_BASIS

#endif /* MINT_FEM_BASIS_HPP_ */
