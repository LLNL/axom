==================
Input Deck Options
==================
--------------
thermal_solver
--------------
----
mesh
----
--
u0
--
-----
kappa
-----
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
   * - filename
     - file for thermal solver
     - Not specified
     - Not specified
     - True
   * - serial
     - serial value
     - Not specified
     - Not specified
     - Not specified
   * - parallel
     - parallel value
     - Not specified
     - Not specified
     - True
   * - order
     - thermal solver order
     - Not specified
     - Not specified
     - True
   * - timestepper
     - thermal solver timestepper
     - Not specified
     - Not specified
     - True
   * - type
     - description for u0 type
     - Not specified
     - Not specified
     - True
   * - func
     - description for u0 func
     - Not specified
     - Not specified
     - True
   * - type
     - description for kappa type
     - Not specified
     - Not specified
     - True
   * - constant
     - description for kappa constant
     - Not specified
     - Not specified
     - True
   * - rel_tol
     - description for solver rel tol
     - Not specified
     - Not specified
     - True
   * - abs_tol
     - description for solver abs tol
     - Not specified
     - Not specified
     - True
   * - print_level
     - description for solver print level
     - Not specified
     - Not specified
     - True
   * - max_iter
     - description for solver max iter
     - Not specified
     - Not specified
     - True
   * - dt
     - description for solver dt
     - Not specified
     - Not specified
     - True
   * - steps
     - description for solver steps
     - Not specified
     - Not specified
     - True
