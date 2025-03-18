// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "MIRApplication.hpp"
#include "axom/mir/utilities/blueprint_utilities.hpp"

#include <conduit_relay_mpi_io_blueprint.hpp>
#include <mpi.h>

namespace bputils = axom::mir::utilities::blueprint;

/*!
 * \brief Create a derived application class that overrides some behaviors for parallel.
 */
class MIRApplicationMPI : public MIRApplication
{
public:
  MIRApplicationMPI() : MIRApplication()
  {
  }
protected:
  /*!
   * \brief Adjust the mesh for parallel by changing the coordinates and
   *        adding state information.
   *
   * \param n_mesh The mesh to modify.
   */
  virtual void adjustMesh(conduit::Node &n_mesh) override
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    n_mesh["state/domain_id"] = rank;

    constexpr int nDomainsPerRow = 8;
    const int domI = rank % nDomainsPerRow;
    const int domJ = rank / nDomainsPerRow;

    // We'll translate the X,Y coordinates for the domain.
    const float xShift = static_cast<float>(domI * gridSize);
    const float yShift = static_cast<float>(domJ * gridSize);
    auto xcView = bputils::make_array_view<float>(n_mesh["coordsets/coords/values/x"]);
    auto ycView = bputils::make_array_view<float>(n_mesh["coordsets/coords/values/y"]);
    for(axom::IndexType i = 0; i < xcView.size(); i++)
    {
      xcView[i] += xShift;
      ycView[i] += yShift;
    }
  }

  /*!
   * \brief Save the mesh to a file using parallel Blueprint routines.
   *
   * \param path The filepath where the file will be saved.
   * \param n_mesh The mesh to be saved.
   */
  virtual void saveMesh(const conduit::Node &n_mesh, const std::string &path) override
  {
#if defined(CONDUIT_RELAY_IO_HDF5_ENABLED)
    std::string protocol("hdf5");
#else
    std::string protocol("yaml");
#endif
    conduit::relay::mpi::io::blueprint::save_mesh(n_mesh, path, protocol, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
};

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
 
  MIRApplicationMPI app;
  int retval = app.initialize(argc, argv);
  if(retval == 0)
  {
    retval = app.execute();
  }

  MPI_Finalize();

  return retval;
}
