.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)


Zero to Axom: Quick install of Axom and Third Party Dependencies
================================================================

The quickest path to install Axom and its dependencies is via `uberenv <https://uberenv.readthedocs.io/en/latest/>`_, a script included in Axom's repo:


.. when on github   git clone --recursive https://github.com/llnl/axom.git?

.. code:: bash

    $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git
    $ cd axom
    $ python scripts/uberenv/uberenv.py --install --prefix="build"


After this completes, ``build/axom-install`` will contain an Axom install.


.. _using_in_your_project:

Using Axom in Your Project
--------------------------

The install includes examples that demonstrate how to use Axom
in a CMake-based and Makefile-based build systems.

CMake-based build system example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. literalinclude:: ../../../examples/using-with-cmake/CMakeLists.txt
   :lines: 35-49

See:  ``examples/axom/using-with-cmake``


Makefile-based build system example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../examples/using-with-make/Makefile
   :lines: 29-35

See: ``examples/axom/using-with-make``
