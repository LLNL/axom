.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

==========================
Iterative Equi-Z Algorithm
==========================

The iterative Equi-Z algorithm provided here is an implementation of the algorithm described 
in the 2010 paper "Visualization and Analysis-Oriented Reconstruction of Material Interfaces" 
by Meredith and Childs.

Functionality for this algorithm is provided by ``computeReconstructedInterfaceIterative()``.
It takes in references to two meshes (the mesh to be processed and the mesh that 
will contain the results of the material interface reconstruction), an integer denoting the 
number of iterations to run for, and a double denoting the percentage difference to modify
the volume fractions by at each iteration.

This algorithm performs reconstruction using the Equi-Z algorithm, but after finishing recontruction,
calculates a percent difference between the resulting element volume fractions and the desired element volume fractions.
It modifies the original element volume fractions by these values and does another reconstruction with 
the Equi-Z algorithm. Since convergence to an optimal value is not guaranteed, it does this for a set number
of iterations that the user specifies.