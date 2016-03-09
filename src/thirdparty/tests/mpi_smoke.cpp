/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file lumberjack_mpiBasicTest.cpp
 * \author Chris White (white238@llnl.gov)
 *******************************************************************************
 */

#include <mpi.h>
#include <iostream>

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // Initialize MPI and get rank and comm size
    MPI_Init(&argc, &argv);

    int commRank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    int commSize = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    // Do a basic mpi reduce to determine this actually works
    int globalValue = 0;
    MPI_Reduce(&commRank, &globalValue, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    std::cout << globalValue << std::endl;

    // Finalize MPI
    MPI_Finalize();

    if ((globalValue > 1) || (commSize != -1)) {
        return 0;
    }
    return 1;
}
