Input
=====

The input to Shroud is a YAML formatted file.
YAML is a human friendly data serialization standard. [yaml]_
Structure is shown through indentation (one or more spaces).  Sequence
items are denoted by a dash, and key value pairs within a map are
separated by a colon::

    library: Tutorial

    types:
      TypeID:
        typedef  : int
        cpp_type : TypeID
    
    functions:
    - decl: void Function1

    classes:
    - name: Class1
      methods:
      - decl: void Method1()

Shroud use curly braces for format strings.
If a string starts with a curly brace YAML
will interpret it as a map/dictionary instead of as part of the
string. To avoid this behavior, strings which start with a curly brace
should be quoted::

    name : "{fmt}"

Strings may be split across several lines by indenting the continued line::

    - decl: void Sum(int len, int *values+dimension+intent(in),
                     int *result+intent(out))

Some values consist of blocks of code.  YAML provides a syntax for 
add multiple lines while preserving newlines::

    C_invalid_name: |
        if (! isNameValid({cpp_var})) {{
            return NULL;
        }}

Note that to insert a literal ``{``, a double brace is required.


Customizing Behavior in the YAML file
-------------------------------------

Fields
^^^^^^

A fields only apply to the type, class or function to which it belongs.
It is not inherited.
For example, *C_name* is a field which is used to explicitly name
a single C wrapper function.  While *C_name_template* is an option which
controls the default value of *C_name*::

    library: testnames

    classes:
      - name: Names
        C_header_filename: foo.h
        C_impl_filename: foo.cpp
        methods:
        -  decl: void method1
           C_name: testmethod1

Annotations
^^^^^^^^^^^

Annotations or attributes apply to specific arguments or results.
They describe semantic behavior for an argument::

    - decl: Class1 *new()  +constructor
    - decl: void delete()  +destructor
    - decl: void Sum(int len, int *values+dimension+intent(in))

Options
^^^^^^^

Options are used to customize the behavior of Shroud.
They are defined in the YAML files as a dictionary.
Options can be defined at the global, class, or function level.
Each level creates a new scope which can access all upper level options.
This allows the user to modifiy behavior for all functions or just a single one::

    options:
      option_a = false
      option_b = false
      option_c = false

    classes:
    - name: class1
      options:
    #    option_a = false     # inherited
         option_b = true
    #    option_c = false     # inherited
      methods:
      - decl: void funtion1
        options:
    #     option_a = false    # inherited
    #     option_b = true     # ihherited
          option_c = true

How code is formatted
---------------------

Format strings contain “replacement fields” surrounded by curly braces
``{}``. Anything that is not contained in braces is considered literal
text, which is copied unchanged to the output. If you need to include
a brace character in the literal text, it can be escaped by doubling:
``{{`` and ``}}``. [Python_Format]_




.. rubric:: Footnotes

.. [Python_Format] https://docs.python.org/2/library/string.html#format-string-syntax

.. [yaml] `yaml.org <http://yaml.org/>`_





