// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file shaping_helpers.hpp
 *
 * \brief Free-standing helper functions in support of shaping query
 */

#ifndef AXOM_QUEST_SHAPING_HELPERS__HPP_
#define AXOM_QUEST_SHAPING_HELPERS__HPP_

#include "mfem.hpp"

namespace axom
{
namespace quest
{
namespace shaping
{
using QFunctionCollection = mfem::NamedFieldsMap<mfem::QuadratureFunction>;
using DenseTensorCollection = mfem::NamedFieldsMap<mfem::DenseTensor>;

/**
 * Utility function to zero out inout quadrature points when a material is replaced
 *
 * When \a shouldReplace is true, we clear all values in \a materialQFunc 
 * that are set in \a shapeQFunc. When it is false, we do the opposite.
 */
void replaceMaterial(mfem::QuadratureFunction* shapeQFunc,
                     mfem::QuadratureFunction* materialQFunc,
                     bool shouldReplace);

/// Utility function to copy inout quadrature point values from one QFunc to another
void copyShapeIntoMaterial(mfem::QuadratureFunction* shapeQFunc,
                           mfem::QuadratureFunction* materialQFunc);

/// Generates a quadrature function corresponding to the mesh positions
void generatePositionsQFunction(mfem::Mesh* mesh,
                                QFunctionCollection& inoutQFuncs,
                                int sampleRes);

/**
 * Compute volume fractions function for shape on a grid of resolution \a gridRes
 * in region defined by bounding box \a queryBounds
 */
void computeVolumeFractions(const std::string& matField,
                            mfem::DataCollection* dc,
                            QFunctionCollection& inoutQFuncs,
                            int outputOrder);

/** 
 * Implements flux-corrected transport (FCT) to convert the inout samples (ones and zeros)
 * to a grid function on the degrees of freedom such that the volume fractions are doubles
 * between 0 and 1 ( \a y_min and \a y_max )
 */
void FCT_project(mfem::DenseMatrix& M,
                 mfem::DenseMatrixInverse& M_inv,
                 mfem::Vector& m,
                 mfem::Vector& x,    // indicators
                 double y_min,       // 0
                 double y_max,       // 1
                 mfem::Vector& xy);  // indicators * rho

namespace unused
{
/// This function is not currently being used, but might be in the near future
void generate_volume_fractions_baseline(mfem::DataCollection* dc,
                                        mfem::QuadratureFunction* inout,
                                        const std::string& name,  // vol_frac
                                        int /*order*/);
}  // end namespace unused

}  // end namespace shaping
}  // end namespace quest
}  // end namespace axom

#endif  // AXOM_QUEST_SHAPING_HELPERS__HPP_
