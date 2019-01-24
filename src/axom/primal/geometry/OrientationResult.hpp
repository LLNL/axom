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

#ifndef PRIMAL_ORIENTATIONRESULT_HPP_
#define PRIMAL_ORIENTATIONRESULT_HPP_

/*!
 * \file
 *
 * \brief Defines the OrientationResult enum which defines possible return
 *  values for orientation tests in between different geometric primitives
 */

namespace axom
{
namespace primal
{

/*!
 * \brief Enumerates possible return values for orientation tests.
 */
enum OrientationResult
{
  ON_BOUNDARY,       /*!< primitive is on the boundary of a primitive      */
  ON_POSITIVE_SIDE,  /*!< primitive is on the positive side of a primitive */
  ON_NEGATIVE_SIDE   /*!< primitive is on the negative side of a primitive */
};

}
}

#endif /* PRIMAL_ORIENTATIONRESULT_HPP_ */
