.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

==================================
Example Output: Input File Options
==================================

--------------
thermal_solver
--------------

.. list-table:: Fields
   :widths: 25 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
     - Value
   * - timestepper
     - thermal solver timestepper
     - quasistatic
     - quasistatic, forwardeuler, backwardeuler
     - |uncheck|
     - quasistatic
   * - order
     - thermal solver order
     - 
     - 1 to 2147483647
     - |check|
     - 2

------
solver
------

Description: This is the solver sub-table in the thermal_solver table

.. list-table:: Fields
   :widths: 25 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
     - Value
   * - dt
     - description for solver dt
     - 1.000000
     - 0.000e+00 to 1.798e+308
     - |check|
     - 1.000000
   * - max_iter
     - description for solver max iter
     - 100
     - 1 to 2147483647
     - |uncheck|
     - 100
   * - print_level
     - description for solver print level
     - 0
     - 0 to 3
     - |check|
     - 0
   * - abs_tol
     - description for solver abs tol
     - 0.000000
     - 0.000e+00 to 1.798e+308
     - |check|
     - 0.000000
   * - steps
     - description for solver steps
     - 1
     - 1 to 2147483647
     - |check|
     - 1
   * - rel_tol
     - description for solver rel tol
     - 0.000001
     - 0.000e+00 to 1.798e+308
     - |uncheck|
     - 0.000001

-----
kappa
-----

.. list-table:: Fields
   :widths: 25 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
     - Value
   * - constant
     - description for kappa constant
     - 
     - 
     - |check|
     - 0.500000
   * - type
     - description for kappa type
     - 
     - constant, function
     - |check|
     - constant

--
u0
--

.. list-table:: Fields
   :widths: 25 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
     - Value
   * - func
     - description for u0 func
     - 
     - 
     - |check|
     - BoundaryTemperature
   * - type
     - description for u0 type
     - constant
     - constant, function
     - |uncheck|
     - function

----
mesh
----

.. list-table:: Fields
   :widths: 25 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
     - Value
   * - parallel
     - 
     - 1
     - 1 to 2147483647
     - |uncheck|
     - 1
   * - serial
     - serial value
     - 1
     - 0 to 2147483647
     - |uncheck|
     - 1
   * - filename
     - file for thermal solver
     - 
     - 
     - |check|
     - /data/star.mesh
