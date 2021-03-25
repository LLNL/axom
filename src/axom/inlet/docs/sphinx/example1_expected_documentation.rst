.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

==================
Example Output: Input File Options
==================

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
   * - timestepper
     - thermal solver timestepper
     - quasistatic
     - quasistatic, forwardeuler, backwardeuler
     - |uncheck|
   * - order
     - thermal solver order
     - 
     - 1 to 2147483647
     - |check|

------
solver
------

Description: This is the solver sub-container in the thermal_solver container

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - dt
     - description for solver dt
     - 1.000000
     - 0.000e+00 to 1.798e+308
     - |check|
   * - max_iter
     - description for solver max iter
     - 100
     - 1 to 2147483647
     - |uncheck|
   * - print_level
     - description for solver print level
     - 0
     - 0 to 3
     - |check|
   * - abs_tol
     - description for solver abs tol
     - 0.000000
     - 0.000e+00 to 1.798e+308
     - |check|
   * - steps
     - description for solver steps
     - 1
     - 1 to 2147483647
     - |check|
   * - rel_tol
     - description for solver rel tol
     - 0.000001
     - 0.000e+00 to 1.798e+308
     - |uncheck|

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
   * - constant
     - description for kappa constant
     - 
     - 
     - |check|
   * - type
     - description for kappa type
     - 
     - constant, function
     - |check|

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
   * - func
     - description for u0 func
     - 
     - 
     - |check|
   * - type
     - description for u0 type
     - constant
     - constant, function
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
   * - parallel
     - 
     - 1
     - 1 to 2147483647
     - |uncheck|
   * - serial
     - serial value
     - 1
     - 0 to 2147483647
     - |uncheck|
   * - filename
     - file for thermal solver
     - 
     - 
     - |check|
