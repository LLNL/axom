Reference
=========

Code
----

The C wrapper code::

    struct s_{C_type_name};
    typedef struct s_{C_type_name} {C_type_name};

    AA_exclass1 * AA_exclass1_new(const char * name);

C implementation::

    {C_return_type} {C_name}({C_prototype})
    AA_exclass1 * AA_exclass1_new(const char * name)
    {
        {C_const}{cpp_class} *{C_this}obj = new {cpp_class}({C_call_list});
        ExClass1 *selfobj = new ExClass1(name);
        return static_cast<AA_exclass1 *>(static_cast<void *>(selfobj));
    }

    void AA_exclass1_delete(AA_exclass1 * self)
    {
        ExClass1 *selfobj = static_cast<ExClass1 *>(
            static_cast<void *>(self));
        delete selfobj;
    }

    int AA_exclass1_increment_count(AA_exclass1 * self, int incr)
    {
        {C_const}{cpp_class} *{C_this}obj =
            static_cast<{C_const}{cpp_class} *>(
                static_cast<{C_const}void *>({C_this}));
        ExClass1 *selfobj = static_cast<ExClass1 *>(
            static_cast<void *>(self));

        {rv_decl} = {CPP_this_call}{method_name}{CPP_template}
            ({C_call_list});
        int rv = selfobj->incrementCount(incr);
        return rv;
    }


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

      interface {F_name_generic}
        module procedure {F_name_impl}
      end interface {F_name_generic}

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

F_name_generic

    method1

F_name_method

    The name of the *F_name_impl* subprogram when used as a
    type procedure.

PY_name_impl

    PY_class1_method1



templates
^^^^^^^^^

Templates are set in options then expanded to assign to the format 
dictionary.

C_name_method_template

    {C_prefix}{lower_class}_{underscore_name}{function_suffix}

C_name_function_template

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

F_name_impl_method_template

    '{lower_class}_{underscore_name}{function_suffix}'

F_name_impl_function_template

    '{underscore_name}{function_suffix}'

F_name_method

    '{underscore_name}{function_suffix}'

F_name_generic

    '{underscore_name}'





C_this

    Name of the C object argument.  Defauls to ``self``.

F_this

    Name of the Fortran object argument.  Defaults to ``obj``.

F_result

    pass

F_derived_member

    pass


Top Level
---------

C_header_filename

   Output file name for header for  wrapper routines.
   Defaults to option *C_header_filename_library_template*.

C_impl_filename

   Output file name for implementation of wrapper routines.
   Defaults to option *C_impl_filename_library_template*.

F_module_name

   Name of Fortran module for this class.
   Defaults to option *F_module_name_library_template*.

F_impl_filename

   Name of Fortran file for functions.
   Defaults to option *F_impl_name_library_template*.


Classes
-------

C_header_filename

   Output file name for header for  wrapper routines.
   Defaults to option *C_header_filename_class_template*.

C_impl_filename

   Output file name for implementation of wrapper routines.
   Defaults to option *C_impl_filename_class_template*.

F_module_name

   Name of Fortran module for this class.
   Defaults to option *F_module_name_class_template*.
   Only used if option *F_module_per_class* is True.

F_impl_filename

   Name of Fortran file for this class.
   Defaults to option *F_impl_name_class_template*.
   Only used if option *F_module_per_class* is True.


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

return_this

   The method returns a reference to ``this``.  This ideom can be used
   to chain calls in C++.  This does not translate to C and Fortran.
   Instead make the return type ``void``.



C_name

    Name of the C wrapper function.
    Defaults to option *C_name_method_template* or
    *C_name_function_template*.

F_C_name

    Name of the Fortran ``BIND(C)`` interface for a C function.
    Defaults to the lower case version of *C_name*.
..    tut_class1_method1

F_name_impl

    Name of the Fortran implementation function.
    Defaults to option *F_name_impl_method_template* or
    *F_name_impl_function_template*.
..    class1_method1



Annotations
-----------

constructor

   Mark method as a constructor.

destructor

   Mark method as a destructor.

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

len_trim

   For a string argument, pass the string address and the result of
   len_trim.
