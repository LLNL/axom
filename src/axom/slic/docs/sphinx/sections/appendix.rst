.. ## Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _sections/appendix:

Appendix
---------

 .. _SlicApplicationCodeExample:

Slic Application Code Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is the complete :ref:`SlicApplicationCodeExample` presented in
the :ref:`sections/getting_started` section. The code can be found in the Axom
source code under ``src/axom/slic/examples/basic/logging.cpp``.

 .. literalinclude:: ../../../examples/basic/logging.cpp
   :start-after: SPHINX_SLIC_BASIC_EXAMPLE_BEGIN
   :end-before: SPHINX_SLIC_BASIC_EXAMPLE_END
   :language: C++
   :linenos:


.. _axomProcessAbort:

axom::utilities::processAbort()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`axomProcessAbort` function gracefully aborts the application by:

#. Calling ``abort()`` if it is a serial application.

#. Calls ``MPI_Abort()`` if the `Axom Toolkit`_ is compiled with MPI and the
   application has initialized MPI, i.e., it's a distributed MPI application.

.. _SlicMacros:

Slic Macros Used in Axom
^^^^^^^^^^^^^^^^^^^^^^^^^
Slic provides a set of convenience macros that can be used to streamline
logging within an application.

.. note::

  The :ref:`SlicMacros` are not the only interface
  to log messages with Slic. They are used within the `Axom Toolkit`_ for
  convenience. Applications or libraries that adopt Slic, typically, use the
  C++ API directly, e.g., call ``slic::logMessage()`` and  wrap the
  functionality in application specific macros to better suit the requirements
  of the application.

.. _CollectiveSlicMacros:

Collective Slic Macros
^^^^^^^^^^^^^^^^^^^^^^^^^
A subset of SLIC macros are collective operations when used with
MPI-aware :ref:`LogStream` instances such as :ref:`SynchronizedStream`
or :ref:`LumberjackStream`.

Additionally, macros such as ``SLIC_WARNING`` and ``SLIC_CHECK`` become collective
operations when certain flags are toggled on or functions are called. Other macros
such as ``SLIC_ERROR`` and ``SLIC_ASSERT`` can be made not collective when certain
functions are called.

The table below details which SLIC macros are collective:

+----------------------------+------------+--------------------------------------------------------------------------+
| Macro                      | Collective | Notes                                                                    |
+============================+============+==========================================================================+
| | ``SLIC_INFO``            | | No       | |                                                                        |
| | ``SLIC_INFO_IF``         | |          | |                                                                        |
| | ``SLIC_INFO_ROOT``       | |          | |                                                                        |
| | ``SLIC_INFO_ROOT_IF``    | |          | |                                                                        |
|                            |            |                                                                          |
| | ``SLIC_ERROR``           | | Yes      | | Not collective after ``slic::disableAbortOnError()`` is called         |
| | ``SLIC_ERROR_IF``        | |          | |                                                                        |
| | ``SLIC_ERROR_ROOT``      | |          | |                                                                        |
| | ``SLIC_ERROR_ROOT_IF``   | |          | |                                                                        |
|                            |            |                                                                          |
| | ``SLIC_WARNING``         | | Yes      | | Collective after ``slic::enableAbortOnWarning()`` is called            |
| | ``SLIC_WARNING_IF``      | |          | |                                                                        |
| | ``SLIC_WARNING_ROOT``    | |          | |                                                                        |
| | ``SLIC_WARNING_ROOT_IF`` | |          | |                                                                        |
|                            |            |                                                                          |
| | ``SLIC_DEBUG``           | | No       | |                                                                        |
| | ``SLIC_DEBUG_IF``        | |          | |                                                                        |
| | ``SLIC_DEBUG_ROOT``      | |          | |                                                                        |
| | ``SLIC_DEBUG_ROOT_IF``   | |          | |                                                                        |
|                            |            |                                                                          |
| | ``SLIC_ASSERT``          | | Yes      | | Not collective after ``slic::disableAbortOnError()`` is called         |
| | ``SLIC_ASSERT_MSG``      | |          | |                                                                        |
|                            |            |                                                                          |
| | ``SLIC_CHECK``           | | Yes      | | Collective after ``slic::debug::checksAreErrors`` is set to true       |
| | ``SLIC_CHECK_MSG``       | |          | |                                                                        |
+----------------------------+------------+--------------------------------------------------------------------------+

Doxygen generated API documentation on Macros can be found here: `SLIC Macros <../../../../doxygen/html/slic__macros_8hpp.html>`_
