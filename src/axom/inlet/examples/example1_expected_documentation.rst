==================
Input Deck Options
==================
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
     - True
   * - serial
     - serial value
     - 
     - 
     - 
   * - parallel
     - 
     - 
     - 
     - True
   * - order
     - thermal solver order
     - 
     - 
     - True
   * - timestepper
     - thermal solver timestepper
     - 
     - 
     - True
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
     - True
   * - func
     - description for u0 func
     - 
     - 
     - True
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
     - True
   * - constant
     - description for kappa constant
     - 
     - 
     - True
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
     - True
   * - abs_tol
     - description for solver abs tol
     - 
     - 
     - True
   * - print_level
     - description for solver print level
     - 
     - 
     - True
   * - max_iter
     - description for solver max iter
     - 
     - 
     - True
   * - dt
     - description for solver dt
     - 
     - 
     - True
   * - steps
     - description for solver steps
     - 
     - 
     - True
