.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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

.. _SLIC_INFO:

SLIC_INFO
"""""""""

The ``SLIC_INFO`` macro logs ``INFO`` messages that consist general
information reported by an application.

The following code snippet illustrates the usage of the ``SLIC_INFO`` macro:

.. code-block:: c++

    SLIC_INFO( "Total number of mesh cells:" <<  N );

SLIC_INFO_IF
""""""""""""

The ``SLIC_INFO_IF`` macro provides the same functionality with the
:ref:`SLIC_INFO` macro, however, it takes one additional argument,
a boolean expression, that allows the application to conditionally log an
``INFO`` message depending on the value of the boolean expression.

For example, the following code snippet illustrates the usage of the
``SLIC_INFO_IF`` macro.

.. code-block:: c++

    SLIC_INFO_IF( (myval >= 0), "[" << myval << "] is positive!" );

In this case, the ``INFO`` message is only logged when the boolean
expression, ``(myval >=0)`` evaluates to ``true``.

.. note::

  The primary intent is to provide a convenience layer and facilitate in a
  cleaner and more compact code style by encapsulating the conditional branching
  logic within a macro.

.. _SLIC_ERROR:

SLIC_ERROR
"""""""""""

The ``SLIC_ERROR`` macro logs ``ERROR`` messages that indicate that the
application has encountered a critical error.

The following code snippet illustrates the basic usage of the
``SLIC_ERROR`` macro:

.. code-block:: c++

    SLIC_ERROR( "jacobian is negative!" );

A stacktrace of the application is appended to all ``ERROR`` messages to
facilitate debugging.

.. note::

  By default, an ``ERROR`` message triggers a call to ``abort()`` the
  application. However, this behavior can be toggled by calling
  ``slic::enableAbortOnError()`` and ``slic::disableAbortOnError()`` accordingly.
  See the `Slic Doxygen API Documentation`_ for more information.

SLIC_ERROR_IF
""""""""""""""

The ``SLIC_ERROR_IF`` provides the same functionality with the :ref:`SLIC_ERROR`
macro, however, it takes one additional argument, a boolean expression, that
allows the application to conditionally log an ``ERROR`` message depending on
the value of the boolean expression.

The following code snippet illustrates the usage of the ``SLIC_ERROR_IF`` macro.

.. code-block:: c++

    SLIC_ERROR_IF( (jacobian < 0.0 + TOL), "jacobian is negative!" );

In this case, the ``ERROR`` message in only logged when the boolean
expression, ``(jacobian < 0.0 + TOL)`` evaluates to ``true``.

.. note::

  The primary intent is to provide a convenience layer and facilitate in a
  cleaner and more compact code style by encapsulating the conditional branching
  logic within a macro.

.. _SLIC_WARNING:

SLIC_WARNING
"""""""""""""

The ``SLIC_WARNING`` macro logs ``WARNING`` messages that indicate that
the application has encountered an error, however, the error is not critical
and the application can proceed.

The following code snippet illustrates the basic usage of the ``SLIC_WARNING``
macro.

.. code-block:: c++

    SLIC_WARNING( "Region [" << ir << "] defined but not used in the problem!" );

SLIC_WARNING_IF
""""""""""""""""

Similarly, the ``SLIC_WARNING_IF`` macro provides the same functionality with
the :ref:`SLIC_WARNING` macro, however, it takes one additional argument,
a boolean expression, that allows the application to conditionally log a
``WARNING`` message depending on the value of the boolean expression.

The following code snippet illustrates the basic usage of the
``SLIC_WARNING_IF`` macro.

.. code-block:: c++

    SLIC_WARNING_IF( (val < 1), "val cannot be less than 1. Forcing value to 1." );
    val = 1;

In this case, the ``WARNING`` message is only logged when the boolean
expression, ``(val < 1)``, evaluates to `` true``.

.. note::

  The primary intent is to provide a convenience layer and facilitate in a
  cleaner and more compact code style by encapsulating the conditional branching
  logic within a macro.

.. _SLIC_DEBUG:

SLIC_DEBUG
"""""""""""

The ``SLIC_DEBUG`` macro logs ``DEBUG`` messages that are intended for
debugging information intended for developers.

The following code snippet illustrates the basic usage of the ``SLIC_DEBUG``
macro

.. code-block:: c++

   SLIC_DEBUG( "Application is running with " << N << " threads!" );

.. warning::

   This macro will log messages only when the `Axom Toolkit`_ is configured and
   built with debug symbols. Consult the `Axom Quick Start Guide`_ for more
   information.

SLIC_DEBUG_IF
""""""""""""""

Similarly, the ``SLIC_DEBUG_IF`` macro provides the same functionality with the
:ref:`SLIC_DEBUG` macro, however, it take one additional argument, a boolean
expression, that allows the application to conditionally log a ``DEBUG``
message depending on the value of the supplied boolean expression.

The following code snippet illustrates the basic usage of the ``SLIC_DEBUG_IF``
macro.

.. code-block:: c++

    SLIC_DEBUG_IF( (value < 0), "value is negative!" );

In this case, the ``DEBUG`` message is only logged when the boolean
expression, ``(value <0)``, evaluates to ``true``.

.. note::

  The primary intent is to provide a convenience layer and facilitate in a
  cleaner and more compact code style by encapsulating the conditional branching
  logic within a macro.

.. warning::

   This macro will log messages only when the `Axom Toolkit`_ is configured and
   built with debug symbols. Consult the `Axom Quick Start Guide`_ for more
   information.

.. _SLIC_ASSERT:

SLIC_ASSERT
"""""""""""

The ``SLIC_ASSERT`` macro is used in a similar manner to the C++ ``assert()``
function call. It evaluates the given expression and logs an ``ERROR``
message if the assertion is not true. The contents of the error message consist
of the supplied expression.

The ``SLIC_ASSERT`` macro is typically used to capture programming errors and
to ensure pre-conditions and post-conditions are satisfied.

The following code snippet illustrates the basic usage of the ``SLIC_ASSERT``
macro.

.. code-block:: c++

    SLIC_ASSERT( data != nullptr );


.. warning::

   This macro will log messages only when the `Axom Toolkit`_ is configured and
   built with debug symbols. Consult the `Axom Quick Start Guide`_ for more
   information.

SLIC_ASSERT_MSG
""""""""""""""""

The ``SLIC_ASSERT_MSG`` macro provides the same functionality with
the :ref:`SLIC_ASSERT` macro, however, it allows the application to supply
a custom message in addition to the boolean expression that is evaluated.

The following code snippet illustrates the basic usage of the ``SLIC_ASSERT_MSG``
macro.

.. code-block:: c++

    SLIC_ASSERT_MSG( data != nullptr, "supplied pointer is null!" );

.. warning::

   This macro will log messages only when the `Axom Toolkit`_ is configured and
   built with debug symbols. Consult the `Axom Quick Start Guide`_ for more
   information.

.. _SLIC_CHECK:

SLIC_CHECK
"""""""""""

The ``SLIC_CHECK`` macro evaluates a given boolean expression, similar to the
:ref:`SLIC_ASSERT` macro. However, in contrast to the :ref:`SLIC_ASSERT` macro,
when the boolean expression evaluates to false, the macro logs a ``WARNING``
instead of an ``ERROR``.

The following code snippet illustrates the basic usage of the ``SLIC_CHECK``
macro.

.. code-block:: c++

    SLIC_CHECK( data != nullptr );


.. warning::

   This macro will log messages only when the `Axom Toolkit`_ is configured and
   built with debug symbols. Consult the `Axom Quick Start Guide`_ for more
   information.

SLIC_CHECK_MSG
""""""""""""""""

The ``SLIC_CHECK_MSG`` macro provides the same functionality with the
:ref:`SLIC_CHECK` macro, however, it allows for the application to supply
a custom message in addition to the boolean expression that is evaluated.

The following code snippet illustrates the basic usage of the ``SLIC_CHECK_MSG``
macro.

.. code-block:: c++

    SLIC_CHECK_MSG( data != nullptr, "supplied pointer is null!" );

.. warning::

   This macro will log messages only when the `Axom Toolkit`_ is configured and
   built with debug symbols. Consult the `Axom Quick Start Guide`_ for more
   information.

.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst
