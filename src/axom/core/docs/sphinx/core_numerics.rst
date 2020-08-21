.. ## Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core numerics
******************************************************

The `axom::numerics` namespace was designed for convenient representation
and use of a mathematical matrix, with accompanying manipulation and 
solver routines.

As an example, the following code shows vector operations.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _vecops_start
   :end-before: _vecops_end
   :language: C++

This example code shows how to construct a matrix.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _matctor_start
   :end-before: _matctor_end
   :language: C++

We can add and multiply matrices, vectors, and scalars, find the determinant,
and extract upper and lower triangular matrices.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _matops_start
   :end-before: _matops_end
   :language: C++

We can use the power method or the Jacobi method to find the eigenvalues and
vectors of a matrix.  Currently, the API of `eigen_solve` is not consistent
with `jacobi_eigensolve` (`eigen_solve` takes a `double` pointer as input
instead of a `Matrix` and the return codes differ); this is an issue we're fixing.
The power method is a stochastic algorithm, computing many matrix-vector
multiplications to produce approximations of a matrix's eigenvalues and vectors.
The Jacobi method is also an iterative algorithm, but it is not stochastic, and
tends to converge much more quickly and stably than other methods.  However,
the Jacobi method is only applicable to symmetric matrices.  In the following
snippet, we show both the power method and the Jacobi method to demonstrate
that they get the same answer.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _eigs_start
   :end-before: _eigs_end
   :language: C++

We can solve a linear system directly or by using LU decomposition and
back-substitution.

.. literalinclude:: ../../examples/core_numerics.cpp
   :start-after: _solve_start
   :end-before: _solve_end
   :language: C++

