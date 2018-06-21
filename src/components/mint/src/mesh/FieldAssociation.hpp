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

#ifndef MINT_FIELDASSOCIATION_HPP_
#define MINT_FIELDASSOCIATION_HPP_

namespace axom
{
namespace mint
{

/*!
 * \brief Enumerates the number of supported associations of field variables
 *  with a corresponding mesh entity, e.g., node, cell-center, face-center, etc.
 */
enum FieldAssociation
{
  ANY_CENTERING = -1,     //!< wild-card used to indicate any centering
  NODE_CENTERED = 0,      //!< used for fields computed at mesh nodes
  CELL_CENTERED,          //!< used for fields computed at cell centers
  FACE_CENTERED,          //!< used for fields computed at face centers
  EDGE_CENTERED,          //!< used for fields computed at edge centers

  NUM_FIELD_ASSOCIATIONS  //!< max number of field associations
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_FIELDASSOCIATION_HPP_ */
