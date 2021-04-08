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
     - 
   * - `kappa`_
     - 
   * - `u0`_
     - 
   * - `mesh`_
     - 
   * - :ref:`timestepper<timestepper1>`
     - thermal solver timestepper
   * - :ref:`order<order1>`
     - polynomial order


.. _timestepper1:

**timestepper**

thermal solver timestepper

  - Default value: quasistatic
  - Valid values: quasistatic, forwardeuler, backwardeuler
  - Optional


.. _order1:

**order**

polynomial order

  - Valid values: 1 to 2147483647
  - Optional



------
solver
------

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - :ref:`dt<dt2>`
     - time step
   * - :ref:`max_iter<max_iter2>`
     - maximum iteration limit
   * - :ref:`print_level<print_level2>`
     - solver print/debug level
   * - :ref:`abs_tol<abs_tol2>`
     - solver absolute tolerance
   * - :ref:`steps<steps2>`
     - number of steps/cycles to take
   * - :ref:`rel_tol<rel_tol2>`
     - solver relative tolerance


.. _dt2:

**dt**

time step

  - Default value: 1.000000
  - Valid values: 0.000e+00 to 1.798e+308
  - Optional


.. _max_iter2:

**max_iter**

maximum iteration limit

  - Default value: 100
  - Valid values: 1 to 2147483647
  - Optional


.. _print_level2:

**print_level**

solver print/debug level

  - Default value: 0
  - Valid values: 0 to 3
  - Optional


.. _abs_tol2:

**abs_tol**

solver absolute tolerance

  - Default value: 0.000000
  - Valid values: 0.000e+00 to 1.798e+308
  - Optional


.. _steps2:

**steps**

number of steps/cycles to take

  - Default value: 1
  - Valid values: 1 to 2147483647
  - Optional


.. _rel_tol2:

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
   * - :ref:`constant<constant3>`
     - thermal conductivity constant
   * - :ref:`type<type3>`
     - description for kappa type


.. _constant3:

**constant**

thermal conductivity constant

  - Optional


.. _type3:

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
   * - :ref:`func<func4>`
     - description for u0 func
   * - :ref:`type<type4>`
     - description for u0 type


.. _func4:

**func**

description for u0 func

  - Optional


.. _type4:

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
   * - :ref:`parallel<parallel5>`
     - 
   * - :ref:`serial<serial5>`
     - number of serial refinements
   * - :ref:`filename<filename5>`
     - mesh filename


.. _parallel5:

**parallel**



  - Default value: 1
  - Valid values: 1 to 2147483647
  - Optional


.. _serial5:

**serial**

number of serial refinements

  - Default value: 1
  - Valid values: 0 to 2147483647
  - Optional


.. _filename5:

**filename**

mesh filename

  - Optional
