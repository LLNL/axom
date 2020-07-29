==================
Input Deck Options
==================
.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

--------------
thermal_solver
--------------

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - order
     - thermal solver order
     - 
     - 1 to 2147483647
     - |check|
   * - timestepper
     - thermal solver timestepper
     - quasistatic
     - quasistatic, forwardeuler, backwardeuler
     - |uncheck|

----
mesh
----

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - filename
     - file for thermal solver
     - 
     - 
     - |check|
   * - serial
     - serial value
     - 1
     - 0 to 2147483647
     - |uncheck|
   * - parallel
     - 
     - 1
     - 1 to 2147483647
     - |uncheck|

--
u0
--

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - type
     - description for u0 type
     - constant
     - constant, function
     - |uncheck|
   * - func
     - description for u0 func
     - 
     - 
     - |check|

-----
kappa
-----

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - type
     - description for kappa type
     - 
     - constant, function
     - |check|
   * - constant
     - description for kappa constant
     - 
     - 
     - |check|

------
solver
------

Description: This is the solver sub-table in the thermal_solver table

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - rel_tol
     - description for solver rel tol
     - 0.000001
     - 0.000e+00 to 1.798e+308
     - |uncheck|
   * - abs_tol
     - description for solver abs tol
     - 0.000000
     - 0.000e+00 to 1.798e+308
     - |check|
   * - print_level
     - description for solver print level
     - 0
     - 0 to 3
     - |check|
   * - max_iter
     - description for solver max iter
     - 100
     - 1 to 2147483647
     - |uncheck|
   * - dt
     - description for solver dt
     - 1.000000
     - 0.000e+00 to 1.798e+308
     - |check|
   * - steps
     - description for solver steps
     - 1
     - 1 to 2147483647
     - |check|
