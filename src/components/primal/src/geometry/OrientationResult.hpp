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

#ifndef PRIMAL_ORIENTEDSIDE_HPP_
#define PRIMAL_ORIENTEDSIDE_HPP_

/*!
 * \file
 *
 * \brief Defines the OrientedSide enum which defines the possible return values
 *  from orientation checks.
 */

namespace axom
{
namespace primal
{

/*!
 * \brief Enumerates possible return values for orientation tests.
 */
enum OrientedSide
{
  ON_BOUNDARY,         /*!< point is on boundary of the given primitive      */
  ON_POSITIVE_SIDE,    /*!< point is on positive side of the given primitive */
  ON_NEGATIVE_SIDE     /*!< point is on negative side of the given primitive */
};

}
}

#endif /* PRIMAL_ORIENTEDSIDE_HPP_ */
