.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

===========================================
MFEM Coefficient Output: Input file Options
===========================================

---
bcs
---


-------------------
Container contents:
-------------------

Description: List of boundary conditions

The input schema defines an array of this table.
For brevity, only one instance is displayed here.

.. list-table:: Functions
   :widths: 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Function Name
     - Description
     - Signature
     - Required
   * - coef
     - The function representing the BC coefficient
     - Double(Vector, Double)
     - |uncheck|
   * - vec_coef
     - The function representing the BC coefficient
     - Vector(Vector, Double)
     - |uncheck|

-----
attrs
-----


-------------------
Container contents:
-------------------

Description: List of boundary attributes

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - 1
     - 
     - 
     - 
     - |uncheck|
   * - 2
     - 
     - 
     - 
     - |uncheck|
   * - 3
     - 
     - 
     - 
     - |uncheck|
