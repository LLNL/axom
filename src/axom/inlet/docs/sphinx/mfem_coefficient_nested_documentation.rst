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
     - 



---
bcs
---




--------------------
Collection contents:
--------------------



List of boundary conditions

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `attrs`_
     - 
   * - :ref:`coef<coef3>`
     - The function representing the BC coefficient
   * - :ref:`vec_coef<vec_coef3>`
     - The function representing the BC coefficient


.. _coef3:

**coef**

The function representing the BC coefficient

  - Signature: Double(Vector, Double)
  - Optional


.. _vec_coef3:

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
