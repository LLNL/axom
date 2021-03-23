.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

===================================
Nested Structs Output (nested mode)
===================================

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `shapes`_
     - 



------
shapes
------




--------------------
Collection contents:
--------------------



.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `geometry`_
     - Geometric information on the shape
   * - `material`_
     - Material of the shape
   * - `name`_
     - Name of the shape


.. _material:

**material**

Material of the shape

  - Valid values: steel, wood, plastic
  - Optional


.. _name:

**name**

Name of the shape




--------
geometry
--------

Geometric information on the shape

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `operators`_
     - List of shape operations to apply
   * - `start_dimensions`_
     - Dimension in which to begin applying operations
   * - `units`_
     - Units for length
   * - `path`_
     - Path to the shape file
   * - `format`_
     - File format for the shape


.. _start_dimensions:

**start_dimensions**

Dimension in which to begin applying operations

  - Default value: 3


.. _units:

**units**

Units for length

  - Default value: cm
  - Valid values: cm, m
  - Optional


.. _path:

**path**

Path to the shape file



.. _format:

**format**

File format for the shape




---------
operators
---------




--------------------
Collection contents:
--------------------

List of shape operations to apply



List of shape operations to apply

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - `slice`_
     - Options for a slice operation
   * - center
     - Center of rotation
   * - `translate`_
     - Translation vector
   * - axis
     - Axis on which to rotate
   * - `rotate`_
     - Degrees of rotation


.. _rotate:

**rotate**

Degrees of rotation

  - Valid values: -1.800e+02 to 1.800e+02
  - Optional



-----
slice
-----

Options for a slice operation

.. list-table::
   :widths: 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Name
     - Description
   * - origin
     - Origin for the slice operation
   * - `z`_
     - z-axis point to slice on
   * - `y`_
     - y-axis point to slice on
   * - `x`_
     - x-axis point to slice on


.. _z:

**z**

z-axis point to slice on



.. _y:

**y**

y-axis point to slice on



.. _x:

**x**

x-axis point to slice on




---------
translate
---------




--------------------
Collection contents:
--------------------

Translation vector
