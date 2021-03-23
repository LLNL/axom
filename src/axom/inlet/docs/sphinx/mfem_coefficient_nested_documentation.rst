.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

=====================================
MFEM Coefficient Output (nested mode)
=====================================
.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `bcs`_
     - List of boundary conditions



---
bcs
---




--------------------
Collection contents:
--------------------

List of boundary conditions



List of boundary conditions

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `attrs`_
     - List of boundary attributes
   * - `coef`_
     - The function representing the BC coefficient
   * - `vec_coef`_
     - The function representing the BC coefficient


.. _coef:

**coef**

The function representing the BC coefficient

  - Signature: Double(Vector, Double)
  - Optional


.. _vec_coef:

**vec_coef**

The function representing the BC coefficient

  - Signature: Vector(Vector, Double)
  - Optional



-----
attrs
-----




--------------------
Collection contents:
--------------------

List of boundary attributes
