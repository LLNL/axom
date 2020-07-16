==================
Input Deck Options
==================
.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

--------------
thermal_solver
--------------


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
     - Range
     - Required
   * - filename
     - file for thermal solver
     - 
     - 
     - |check|
   * - serial
     - serial value
     - 
     - 
     - |uncheck|
   * - parallel
     - 
     - 
     - 
     - |uncheck|
   * - order
     - thermal solver order
     - 
     - 
     - |check|
   * - timestepper
     - thermal solver timestepper
     - 
     - 
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
     - Range
     - Required
   * - type
     - description for u0 type
     - 
     - 
     - |check|
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
     - Range
     - Required
   * - type
     - description for kappa type
     - 
     - 
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
     - Range
     - Required
   * - rel_tol
     - description for solver rel tol
     - 
     - 
     - |uncheck|
   * - abs_tol
     - description for solver abs tol
     - 
     - 
     - |check|
   * - print_level
     - description for solver print level
     - 
     - 
     - |check|
   * - max_iter
     - description for solver max iter
     - 
     - 
     - |uncheck|
   * - dt
     - description for solver dt
     - 
     - 
     - |check|
   * - steps
     - description for solver steps
     - 
     - 
     - |check|
