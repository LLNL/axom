// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PSTLREADER_HPP_
#define QUEST_PSTLREADER_HPP_

#include "axom/core/Macros.hpp"
#include "axom/quest/stl/STLReader.hpp"  // base class

#include "mpi.h"

namespace axom
{
namespace quest
{
class PSTLReader : public STLReader
{
public:
  /*!
   * \brief Constructor.
   * \param [in] comm user-supplied MPI communicator.
   */
  PSTLReader(MPI_Comm comm);

  /*!
   * \brief Destructor.
   */
  virtual ~PSTLReader();

  /*!
   * \brief Reads in an STL file to all ranks in the associated communicator.
   * \note Rank 0 reads in the STL mesh file and broadcasts the data the rest
   *  of the ranks.
   * \return status set to zero on success; set to a non-zero value otherwise.
   */
  virtual int read() final override;

private:
  /*!
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent its use in application code.
   */
  PSTLReader() : m_comm(MPI_COMM_NULL), m_my_rank(0) {};

  MPI_Comm m_comm; /*!< MPI communicator */
  int m_my_rank;   /*!< MPI rank ID      */

  DISABLE_COPY_AND_ASSIGNMENT(PSTLReader);
  DISABLE_MOVE_AND_ASSIGNMENT(PSTLReader);
};

}  // end namespace quest
}  // end namespace axom

#endif /* QUEST_PSTLREADER_HPP_ */
