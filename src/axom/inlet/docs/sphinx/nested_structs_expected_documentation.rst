.. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX
.. |check|      unicode:: U+2611 .. CHECKED BOX

=========================================
Nested Structs Output: Input file Options
=========================================

------
shapes
------


--------------------
Collection contents:
--------------------

The input schema defines a collection of this container.
For brevity, only one instance is displayed here.

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - material
     - Material of the shape
     - 
     - steel, wood, plastic
     - |check|
   * - name
     - Name of the shape
     - 
     - 
     - |check|

--------
geometry
--------

Description: Geometric information on the shape

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - start_dimensions
     - Dimension in which to begin applying operations
     - 3
     - 
     - |uncheck|
   * - units
     - Units for length
     - cm
     - cm, m
     - |uncheck|
   * - path
     - Path to the shape file
     - 
     - 
     - |check|
   * - format
     - File format for the shape
     - 
     - 
     - |check|

---------
operators
---------


--------------------
Collection contents:
--------------------

Description: List of shape operations to apply

The input schema defines a collection of this container.
For brevity, only one instance is displayed here.

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - rotate
     - Degrees of rotation
     - 
     - -1.800e+02 to 1.800e+02
     - |uncheck|

-----
slice
-----

Description: Options for a slice operation

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - z
     - z-axis point to slice on
     - 
     - 
     - |uncheck|
   * - y
     - y-axis point to slice on
     - 
     - 
     - |uncheck|
   * - x
     - x-axis point to slice on
     - 
     - 
     - |uncheck|

---------
translate
---------


--------------------
Collection contents:
--------------------

Description: Translation vector

.. list-table:: Fields
   :widths: 25 25 25 25 25
   :header-rows: 1
   :stub-columns: 1

   * - Field Name
     - Description
     - Default Value
     - Range/Valid Values
     - Required
   * - 0
     - 
     - 
     - 
     - |uncheck|
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
