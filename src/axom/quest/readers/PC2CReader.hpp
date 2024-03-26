// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PC2CREADER_HPP_
#define QUEST_PC2CREADER_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_C2C
  #error PC2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/core/Macros.hpp"
#include "axom/quest/readers/C2CReader.hpp"  // base class

#include "mpi.h"

namespace axom
{
namespace quest
{
class PC2CReader : public C2CReader
{
public:
  PC2CReader() = delete;
  PC2CReader(MPI_Comm comm);

  virtual ~PC2CReader() = default;

  /*!
   * \brief Reads in a C2C file to all ranks in the associated communicator
   *
   * \note Rank 0 reads in the C2C mesh file and broadcasts to the other ranks
   * \return status set to zero on success; non-zero otherwise
   */
  int read() final override;

private:
  /// Utility function to broadcast a vector of primitive types
  template <typename T>
  void bcastVector(std::vector<T>& vec);

private:
  MPI_Comm m_comm {MPI_COMM_NULL};
  int m_my_rank {0};

  DISABLE_COPY_AND_ASSIGNMENT(PC2CReader);
  DISABLE_MOVE_AND_ASSIGNMENT(PC2CReader);
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_PC2CREADER_HPP_
