Specifying Shapes
=================

"Shaping", or "painting", is the process of adding non-conformal material
regions to a mesh. Traditionally, this has been done in code-specific formats
by each code that provides such a capability. Axom's Klee library provides
a way to read shape specifications in YAML files and apply the specified
geometry to a mesh.

Basics
------
Shapes in Klee are specified in YAML. A basic file consists of a list of
shapes. Each one specifies its :code:`name` and :code:`material`,
as well as a description of its geometry.

In addition to the shapes themselves, a file must specify the number of
dimensions of the shapes (only 2 and 3 are allowed). This will be important
later when specifying transformations and slices.

.. code-block:: yaml

    dimensions: 3

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: stl
          path: wheel.stl
      - name: windshield
        material: glass
        geometry:
          format: stl
          path: windshield.stl


The above example describes a series of 3D shapes. The first is a wheel
made of steel. Its geometry is specified in an STL file named :code:`wheel.stl`.
The second shape is named "windshield", is made of "glass", and its geometry
is specified in :code:`windshield.stl`. Note that Klee does not specify
what a particular material means. A material is simply a label which can
be used by a host code to apply properties to part of a mesh.

Paths
*****
The paths specified in shapes are specified either as absolute paths
(begin with a :code:`/`), or as relative paths. Absolute paths are evaluated
as absolute paths in the file systems. Relative paths are evaluated relative
to the YAML file (not the current working directory). For example, in the
file :code:`/path/to/my_shapes.yaml`, the table below illustrates how
different paths would be specified.

+---------------------------------+-----------------------------+
| Path in /path/to/my_shapes.yaml | Path in file system         |
+=================================+=============================+
| /an/absolute/path.stl           | /an/absolute/path.stl       |
+---------------------------------+-----------------------------+
| just_a_file.stl                 | /path/to/just_a_file.stl    |
+---------------------------------+-----------------------------+
| dir/and/a/file.stl              | /path/to/dir/and/a/file.stl |
+---------------------------------+-----------------------------+

Changing Dimensions on a Per-Shape Basis
****************************************
Sometimes it is useful to bring in a geometry file in a different
dimensionality than the one you are working in. For example, you may be
working in 2D, but may want to bring in a 3D file and then slice it
(slices are described below). To do this, you need to specify the
:code:`initial_dimensions` field on the :code:`geometry` of a :code:`shape`.

.. code-block:: yaml

    dimensions: 2

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: stl
          path: wheel.stl
          initial_dimensions: 3
          operators:
            - slice:
                x: 10

Overlay Rules
-------------
Shapes are added to meshes in the order in which they appear in the YAML
file. By default, each one replaces all materials that occupy the space
specified by its geometry file. This can be overridden by using the
:code:`replaces` and :code:`does_not_replace` properties.

.. code-block:: yaml

    dimensions: 3

    shapes:
      - name: wheel
        material: steel
        replaces: [rubber, air]
        geometry:
          format: stl
          path: wheel.stl
      - name: windshield
        does_not_replace: [steel]
        material: glass
        geometry:
          format: stl
          path: windshield.stl

In the example above, the wheel would only replace rubber and air. Any other
materials that happen to be in the same space as it would be left untouched.
The windshield would replace everything except steel. It is an error to
specify both :code:`replaces` and :code:`does_not_replace`.

Operators
---------

When assembling complex geometries, it is often the case that different parts
are specified in different coordinate systems. For example, a description
of the wheel of a car might be specified around its center, not its position
relative to the rest of the car. To help with this, Klee provides a mechanism
to apply transformations to shapes.

.. code-block:: yaml

    dimensions: 3

    shapes:
      - name: windshield
        material: glass
        geometry:
          format: stl
          path: windshield.stl
          operators:
            - rotate: 90
              axis: [0, 1, 0]
              center: [0, 0, -10]
            - translate: [10, 20, 30]

In the example above, the wheel is rotated 90 degrees counterclockwise
around an axis centered at the point :code:`(0, 0, -10)` and pointing in the
direction of the vector :code:`(0, 1, 0)`. It is then translated by the
vector :code:`(10, 20, 30)`.

Supported Operators
*******************
The supported operators are listed below. Unless otherwise specified,
the only difference between the 2D and 3D versions are that whenever points
or vectors are expected, the points and vectors must be of the dimensionality
specified by the shape file.

Operators take the form of :code:`operator_name: value`, where
:code:`operator_name` is the name of the operator, and
:code:`value` is the value specifying the parameters of the operation.
Operators may also have additional required or optional parameters.

* Translations

  :description: Translate the shape by a given vector.
  :name: :code:`translate`
  :value: a vector specifying the amount by which to translate the shape
  :example:
    ::

        # Translate by vector (1, 2, 3)
        translate: [1, 2, 3]

* Rotations

  :description: Rotate the shape by a given amount around a specified axis
  :name: :code:`rotate`
  :value: an angle, in degrees by which the shape will be rotated
    counterclockwise.
  :additional required parameters:
    :axis: (3D only) the axis of rotation
  :optional arguments:
    :center: a point specifying the center of rotation
  :example:
    ::

        # Rotate 45 degrees counterclockwise around the ray passing through
        # the point (1, 2, 3) and pointing in the direction of the vector
        # (4, 5, 6)
        rotate: 45
        center: [1, 2, 3]
        axis: [4, 5, 6]

* Scaling

  :description: Scale the shape by a specified amount
  :name: :code:`scale`
  :value: a vector specifying the amount by which to scale in each dimension,
    or a single value specifying by which to scale in all dimensions
  :example:
    ::

        # Scale by 2x in the x direction 0.5x in y, and 1.5x in z
        scale: [2.0, 0.5, 1.5]

* Arbitrary Affine Matrices

  :description: Apply an arbitrary affine transformation to a shape
  :name: :code:`affine`
  :value: a vector containing either 6 (for 2D) or 12 (for 3D) values. In 2D,
   the vector :code:`(a, b, c, d, e, f)` maps to the matrix
   :code:`((a, b, c), (d, e, f), (0, 0, 1))`. In 3D, the vector
   :code:`(a, b, c, d, e, f, g, h, i, j, k, l)` maps to the matrix
   :code:`((a, b, c, d), (e, f, g, h), (i, j, k, l), (0, 0, 0, 1))`.
  :example:
    ::

        # Apply the matrix ((1, 2, 3, 4,), (5, 6, 7, 8), (9, 10, 11, 12))
        affine: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

* Slices

  :description: Slice a 3D object and convert it into a 2D object. This is
    accomplished by defining a cut plane which will be used to determine
    what slice of the geometry to take. In addition, a point on the plane
    is picked as the new origin, and a vector is used to specify how the
    plane should be oriented with the 2D axes.
  :name: :code:`slice`
  :value: an object with the the following properties

    :origin: the point to use as the origin of the new coordinate system
    :normal: a vector normal to the slice plane
    :up: a vector which will be mapped to the positive Y direction.
  :optional arguments:
    :x: a single value specifying that the cut plane perpendicular to the
      x-axis at this value. See defaults table below.
    :y: a single value specifying that the cut plane perpendicular to the
      y-axis at this value. See defaults table below.
    :z: a single value specifying that the cut plane perpendicular to the
      z-axis at this value. See defaults table below.

    If a plane is specified by just giving "x", "y", or "z", then the origin,
    normal, and up vectors are given the default values specified
    in the table below. They can be overridden so long as the origin is still
    on the plane, and the normal is a multiple of the default normal.

    +------------------+-----------------------+-------------------+-------------------+
    | Usage            | :code:`origin`        | :code:`normal`    | :code:`up`        |
    +==================+=======================+===================+===================+
    | :code:`x: <val>` | :code:`(<val>, 0, 0)` | :code:`(1, 0, 0)` | :code:`(0, 0, 1)` |
    +------------------+-----------------------+-------------------+-------------------+
    | :code:`y: <val>` | :code:`(0, <val>, 0)` | :code:`(0, 1, 0)` | :code:`(1, 0, 0)` |
    +------------------+-----------------------+-------------------+-------------------+
    | :code:`z: <val>` | :code:`(0, 0, <val>)` | :code:`(0, 0, 1)` | :code:`(0, 1, 0)` |
    +------------------+-----------------------+-------------------+-------------------+

  :example:
    ::

        # Cut a 3D object with a plane that passes through the point
        # [10, 20, 30] and is normal to the vector [4, 5, 6]. The vector
        # [-5, 4, 0] will be mapped to the positive Y axis. [10, 20, 30] will
        # be mapped to the origin.
        slice:
          origin: [10, 20, 30]
          normal: [4, 5, 6]
          up: [-5, 4, 0]

Named Operators
***************
It can often be useful to name and reuse operators. For example, you may
have several parts of an assembly specified in one coordinate system
that you then need to transform to another. To enable reuse, we provide
support for named operators.

Named operators are specified via the top-level :code:`named_operators`
object. This is a list where each entry has the following values:

:name (required): the name of the operator. This is how it is referenced later.
:value (required): A list of operators. This is identical to the
  :code:`operators` entry in the :code:`geometry` object of a :code:`shape`.
:initial_dimensions (optional): the number of initial dimensions of the
  operator. Must be 2 or 3. If not specified, the number of dimensions of the
  document is used.

The example below demonstrates how to create and then use a named operator.
Note that :code:`ref` is just one entry in the list of operators, and
additional operators may be used.

.. code-block:: yaml

    dimensions: 2

    shapes:
      - name: wheel
        material: steel
        geometry:
          format: stl
          path: wheel.stl
          operators:
            - ref MyFirstOperator
            - rotate: 90

    named_operators:
      - name: MyFirstOperator
        value:
          - translate: [10, 20]
          - scale: 1.5


In addition to using :code:`ref` in an individual shape's operators, you
can also use it in other named operators. The only restriction is that it
be defined in the list before it is used.

.. code-block:: yaml

    dimensions: 2

    named_operators:
      - name: SomeOperator
        value:
          - translate: [10, 20]
          - scale: 1.5
      - name: AnotherOperator
        value:
          - rotate: 90
          - ref: SomeOperator
