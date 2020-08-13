
Lumberjack User Documentation
=============================

Lumberjack, named because it cuts down logs, is a C++ library that
provides scalable logging while reducing the amount of messages
written out the screen or file system.


.. raw:: html

    <h3>Introduction</h3>

Lumberjack was created to provide scalable logging with a simple programming
model while allowing developers to customize its behavior. It is named Lumberjack
because it cuts down logs. It uses MPI and a scalable binary tree reduction
scheme to combine duplicate messages and limit output to only the root node.


.. raw:: html

    <h3>Requirements</h3>

* MPI - MPI is fundamental to Lumberjack and without MPI, Lumberjack is not useful.


.. raw:: html

    <h3>Code Guarding</h3>

You tell if Axom was built with Lumberjack enabled by using the following
include and compiler define:

.. code-block:: c

    #include "axom/config.hpp"
    #ifdef AXOM_USE_LUMBERJACK
        // Lumberjack code
    #endif


.. raw:: html

    <h3>Contents</h3>


.. toctree::
   :maxdepth: 2

   quick_start
   core_concepts
   lumberjack_classes

.. raw:: html

    <h3>Additional links</h3>


* `API documentation <../../../../doxygen/html/lumberjacktop.html>`_
* `Axom main docs <../../../../index.html>`_
