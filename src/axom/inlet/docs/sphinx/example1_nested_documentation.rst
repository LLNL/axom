.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

================================================
Example Output: Input File Options (nested mode)
================================================

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `thermal_solver`_
     - 



--------------
thermal_solver
--------------

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `solver`_
     - linear equation solver options
   * - `kappa`_
     - 
   * - `u0`_
     - 
   * - `mesh`_
     - 
   * - `timestepper`_
     - thermal solver timestepper
   * - `order`_
     - polynomial order


.. _timestepper:

**timestepper**

thermal solver timestepper

  - Default value: quasistatic
  - Valid values: quasistatic, forwardeuler, backwardeuler
  - Optional


.. _order:

**order**

polynomial order

  - Valid values: 1 to 2147483647
  - Optional



------
solver
------

linear equation solver options

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `dt`_
     - time step
   * - `max_iter`_
     - maximum iteration limit
   * - `print_level`_
     - solver print/debug level
   * - `abs_tol`_
     - solver absolute tolerance
   * - `steps`_
     - number of steps/cycles to take
   * - `rel_tol`_
     - solver relative tolerance


.. _dt:

**dt**

time step

  - Default value: 1.000000
  - Valid values: 0.000e+00 to 1.798e+308
  - Optional


.. _max_iter:

**max_iter**

maximum iteration limit

  - Default value: 100
  - Valid values: 1 to 2147483647
  - Optional


.. _print_level:

**print_level**

solver print/debug level

  - Default value: 0
  - Valid values: 0 to 3
  - Optional


.. _abs_tol:

**abs_tol**

solver absolute tolerance

  - Default value: 0.000000
  - Valid values: 0.000e+00 to 1.798e+308
  - Optional


.. _steps:

**steps**

number of steps/cycles to take

  - Default value: 1
  - Valid values: 1 to 2147483647
  - Optional


.. _rel_tol:

**rel_tol**

solver relative tolerance

  - Default value: 0.000001
  - Valid values: 0.000e+00 to 1.798e+308
  - Optional



-----
kappa
-----

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `constant`_
     - thermal conductivity constant
   * - `type`_
     - description for kappa type


.. _constant:

**constant**

thermal conductivity constant



.. _type:

**type**

description for kappa type

  - Valid values: constant, function
  - Optional



--
u0
--

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `func`_
     - description for u0 func
   * - `type`_
     - description for u0 type


.. _func:

**func**

description for u0 func



.. _type:

**type**

description for u0 type

  - Default value: constant
  - Valid values: constant, function
  - Optional



----
mesh
----

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `parallel`_
     - 
   * - `serial`_
     - number of serial refinements
   * - `filename`_
     - mesh filename


.. _parallel:

**parallel**



  - Default value: 1
  - Valid values: 1 to 2147483647
  - Optional


.. _serial:

**serial**

number of serial refinements

  - Default value: 1
  - Valid values: 0 to 2147483647
  - Optional


.. _filename:

**filename**

mesh filename
