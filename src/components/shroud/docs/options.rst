Reference
=========

Code
----

The template for Fortran code showing names which may 
be controlled directly by the input file::

    module {F_module_name}

      type {F_derived_name}
        type(C_PTR) {F_derived_member}
      contains
        procedure :: {F_name_method} => {F_name_impl}
        generic :: {F_name_generic} => {F_name_method}, ...
      end type {F_derived_name}

      interface
        subroutine {F_C_name} bind(C, name="{C_name}")
          ...
        end subroutine {F_C_name}
      end interface

    contains

      subroutine {F_name_impl}
        ...
      end subroutine {F_name_impl}

    end module {F_module_name}


Options
-------

Options control the behavior of shroud.  They are nested so that 
options applied at the outer most level apply to inner groups.
Options set on a class apply to all methods in that class.


library

  The name of the library.
  Used to name output files and modules.
  The first three letters are used as the default for **C_prefix**.

C_prefix

  Prefix added to name of generated C routines.
  The prefix helps to ensure unique global names.

cpp_header

  C++ header file name.

F_string_result_as_arg

  Function with return a ``char *`` will instead by converted to
  subroutine which require an additional argument for the result.

F_string_len_trim

  For each function with a ``std::string`` argument, create another C
  function which accepts a buffer and length.  The C wrapper will call
  the ``std::string`` constructor, instead of the Fortran wrapper
  creating a ``NULL`` terminated string using ``trim``.  This avoids
  copying the string in the Fortran wrapper.

namespace

  Blank delimited list of namespaces for **cpp_header**.




Names
-----

This section describes options used to name generated functions and
methods.

Each method maintains a dictionary of names which can be
used as part of computed names.  This is refered to as the
format dictionary.

method_name

    The C++ name of the function is extracted from the ``decl`` field.
    This name is typically camel case.

underscore_name

    method_name converted from camel case into snake case.
    ``getName`` becomes ``get_name``.

function_suffix

    Function field 'function_suffix'.

C_name

    Name of C wrapper function.

F_C_name

    tut_class1_method1", 

F_name_generic

    method1

F_name_impl

    class1_method1

F_name_method

    method1

PY_name_impl

    PY_class1_method1



templates
^^^^^^^^^

Templates are set in options then expanded to assign to the format 
dictionary.

C_name_method_template

    {C_prefix}{lower_class}_{underscore_name}{function_suffix}

C_name_function_tempate

    {C_prefix}{underscore_name}{function_suffix}



F_C_name

    defaults to C_name.lower() - tut_class1_method1

F_name_generic_template

    defaults to '{underscore_name'} - method1

F_name_impl

    class - {lower_class}_{underscore_name}{function_suffix}
    function - '{underscore_name}{function_suffix}'
    class1_method1

F_name_method

    {underscore_name}{function_suffix}
    method1

PY_name_impl

    PY_class1_method1




C_header_filename_library_template

   'wrap{library}.h'

C_impl_filename_library_template

    'wrap{library}.cpp'

C_header_filename_class_template

    'wrap{cpp_class}.h'

C_impl_filename_class_template

    'wrap{cpp_class}.cpp'


F_module_name_library_template

    '{lower_library}_mod'

F_impl_filename_library_template

    'wrapf{lower_library}.f'

F_module_name_class_template

    '{lower_class}_mod'

F_impl_filename_class_template

    'wrapf{cpp_class}.f'








C_this

    Name of the C object argument.  Defauls to ``self``.

F_this

    Name of the Fortran object argument.  Defaults to ``obj``.

F_result

    pass

F_derived_member

    pass









Functions
---------

Each function can define fields to define the function
and how it should be wrapped.  These fields apply only
to a single function i.e. they are not inherited.


decl

   Function declaration.
   Parsed to extract function name, type and arguments descriptions.

function_suffix

   Suffix to append to the end of generated name.


Annotations
-----------

constructor

   Mark function as a constructor.

pure

   Sets the Fortran PURE attribute.

dimension

   Sets the Fortran DIMENSION attribute.
   Pointer argument should be passed through since it is an
   array.  *value* must be *False*
   If set without a value, it defaults to ``(*)``.

value

   If true, pass-by-value; else, pass-by-reference.

intent

   Valid valid values are ``in``, ``out``, ``inout``.
   If the argument is ``const``, the default is ``in``.

ptr

   Argument is a pointer

reference

   Argument is a reference

default

   If set the ``optional`` keyword is added to the Fortran interface.
