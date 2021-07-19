// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/readers/PC2CReader.hpp"

#ifndef AXOM_USE_C2C
  #error PC2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/slic.hpp"

namespace axom
{
namespace quest
{
namespace
{
constexpr int READER_SUCCESS = 0;
constexpr int READER_FAILED = -1;

}  // end anonymous namespace

//------------------------------------------------------------------------------
PC2CReader::PC2CReader(MPI_Comm comm) : m_comm(comm)
{
  MPI_Comm_rank(m_comm, &m_my_rank);
}

//------------------------------------------------------------------------------
int PC2CReader::read()
{
  SLIC_ASSERT(m_comm != MPI_COMM_NULL);

  // Clear internal data-structures
  this->clear();

  int rc = -1;  // return code

  switch(m_my_rank)
  {
  case 0:

    rc = C2CReader::read();
    if(rc == READER_SUCCESS)
    {
      int numNURBS = m_nurbsData.size();
      MPI_Bcast(&numNURBS, 1, axom::mpi_traits<int>::type, 0, m_comm);

      for(int i = 0; i < numNURBS; ++i)
      {
        auto& nd = m_nurbsData[i];

        // send the order
        MPI_Bcast(&nd.order, 1, axom::mpi_traits<unsigned>::type, 0, m_comm);

        // send the knots, weights and spline interp params vectors
        bcastVector<double>(nd.knots);
        bcastVector<double>(nd.weights);
        bcastVector<double>(nd.splineInterpolationPointParameters);

        // pack R and Z coordinates of control points into arrays and bcast
        // note: when converting to NURBSData, all units are converted to m_lengthUnit
        // so we do not need to transmit this to the other ranks
        std::vector<double> r, z;
        r.reserve(nd.controlPoints.size());
        z.reserve(nd.controlPoints.size());
        for(const auto& pt : nd.controlPoints)
        {
          r.push_back(pt.getR().getValue());
          z.push_back(pt.getZ().getValue());
        }
        bcastVector<double>(r);
        bcastVector<double>(z);
      }
    }
    else
    {
      MPI_Bcast(&rc, 1, MPI_INT, 0, m_comm);
    }
    break;

  default:

    // Rank 0 broadcasts the number of NURBSData entities, a positive integer, if the
    // C2C file is read successfully, or sends a READER_FAILED flag, indicating
    // that the read was not successful.
    int numNURBS = -1;
    MPI_Bcast(&numNURBS, 1, axom::mpi_traits<int>::type, 0, m_comm);

    if(numNURBS != READER_FAILED)
    {
      rc = READER_SUCCESS;

      m_nurbsData.reserve(numNURBS);

      // Receive and reconstruct each NURBSData
      for(int i = 0; i < numNURBS; ++i)
      {
        c2c::NURBSData nd;

        // receive the NURBS order
        MPI_Bcast(&nd.order, 1, axom::mpi_traits<unsigned>::type, 0, m_comm);

        // receive the knots, weights and spline interp params vectors
        bcastVector<double>(nd.knots);
        bcastVector<double>(nd.weights);
        bcastVector<double>(nd.splineInterpolationPointParameters);

        // get the points data and reconstruct the control points
        std::vector<double> r, z;
        bcastVector<double>(r);
        bcastVector<double>(z);
        int ptSize = r.size();
        for(int j = 0; j < ptSize; ++j)
        {
          const c2c::Length rlen {r[j], m_lengthUnit};
          const c2c::Length zlen {z[j], m_lengthUnit};
          nd.controlPoints.emplace_back(c2c::Point(rlen, zlen));
        }

        m_nurbsData.emplace_back(nd);
      }
    }
  }

  return (rc);
}

template <typename T>
void PC2CReader::bcastVector(std::vector<T>& vec)
{
  int sz = -1;

  switch(m_my_rank)
  {
  case 0:
    sz = vec.size();
    MPI_Bcast(&sz, 1, axom::mpi_traits<int>::type, 0, m_comm);
    MPI_Bcast(vec.data(), sz, axom::mpi_traits<T>::type, 0, m_comm);
    break;
  default:
    MPI_Bcast(&sz, 1, axom::mpi_traits<int>::type, 0, m_comm);
    vec.resize(sz);  // resize to allocate sufficient data
    MPI_Bcast(vec.data(), sz, axom::mpi_traits<T>::type, 0, m_comm);
    break;
  }
}

}  // namespace quest
}  // namespace axom
