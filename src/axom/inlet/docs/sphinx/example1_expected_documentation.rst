.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

==================================
Example Output: Input file Options
==================================

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
     - polynomial order
     - 
     - 1 to 2147483647
     - |check|

------
solver
------

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
     - time step
     - 1.000000
     - 0.000e+00 to 1.798e+308
     - |check|
   * - steps
     - number of steps/cycles to take
     - 1
     - 1 to 2147483647
     - |check|
   * - print_level
     - solver print/debug level
     - 0
     - 0 to 3
     - |check|
   * - rel_tol
     - solver relative tolerance
     - 0.000001
     - 0.000e+00 to 1.798e+308
     - |uncheck|
   * - max_iter
     - maximum iteration limit
     - 100
     - 1 to 2147483647
     - |uncheck|
   * - abs_tol
     - solver absolute tolerance
     - 0.000000
     - 0.000e+00 to 1.798e+308
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
   * - constant
     - thermal conductivity constant
     - 
     - 
     - |check|
   * - type
     - description for kappa type
     - 
     - constant, function
     - |check|

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
   * - filename
     - mesh filename
     - 
     - 
     - |check|
   * - serial
     - number of serial refinements
     - 1
     - 0 to 2147483647
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
