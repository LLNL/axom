.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _sections/appendix:

Appendix
---------

.. _MintApplicationCodeExample:

Mint Application Code Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is the complete :ref:`MintApplicationCodeExample` presented in
the :ref:`sections/quick_introduction` section. The code can be found in the Axom
source code under ``src/axom/mint/examples/mint_uniform_mesh.cpp``.

.. literalinclude:: ../../examples/mint_uniform_mesh.cpp
   :start-after: sphinx_tutorial_basic_example_start
   :end-before: sphinx_tutorial_basic_example_end
   :language: C++
   :linenos:

.. _axomLambdaMacro:

AXOM_LAMBDA Macro
^^^^^^^^^^^^^^^^^

The ``AXOM_LAMBDA`` convenience macro expands to:

 * ``[=]`` capture by value when the `Axom Toolkit`_ is compiled without CUDA.
 * ``[=] __host__ __device__`` when the `Axom Toolkit`_ is compiled with CUDA

.. _rawSidreData:

Raw Sidre Data
^^^^^^^^^^^^^^

.. literalinclude:: sections/raw_sidre_data.txt
   :language: json
   :linenos:
