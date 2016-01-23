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

debug

  Print additional comments in generated files that may 
  be useful for debugging.

library

  The name of the library.
  Used to name output files and modules.
  The first three letters are used as the default for **C_prefix**.

C_prefix

  Prefix added to name of generated C routines.
  The prefix helps to ensure unique global names.

C_proto_type

   XXX  override prototype of generated C function

C_return_type

   XXX   override return type of function

cpp_header

  C++ header file name.

F_string_result_as_arg

  The name of the output argument.
  Function which return a ``char *`` will instead by converted to a
  subroutine which require an additional argument for the result.

F_string_len_trim

  For each function with a ``std::string`` argument, create another C
  function which accepts a buffer and length.  The C wrapper will call
  the ``std::string`` constructor, instead of the Fortran wrapper
  creating a ``NULL`` terminated string using ``trim``.  This avoids
  copying the string in the Fortran wrapper.
  Defaults to *true*.
.. bufferify

F_force_wrapper

  If *true*, always create an explicit Fortran wrapper.
  If *false*, only create the wrapper when there is work for it to do;
  otherwise, call the C function directly.
  For example, a function which only deals with native
  numeric types does not need a wrapper since it can be called
  directly by defining the correct interface.
  The default is *false*.

namespace

  Blank delimited list of namespaces for **cpp_header**.

wrap_c

  If *true*, create C wrappers.
  Defaults to *true*.

wrap_fortran

  If *true*, create Fortran wrappers.
  Defaults to *true*.

wrap_python

  If *true*, create Python wrappers.
  Defaults to *false*.



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

PY_name_impl

    PY_class1_method1



templates
^^^^^^^^^

Templates are set in options then expanded to assign to the format 
dictionary.

C_name_function_template

    {C_prefix}{underscore_name}{function_suffix}

C_name_method_template

    {C_prefix}{lower_class}_{underscore_name}{function_suffix}



F_C_name

    Defaults to C_name.lower() - tut_class1_method1

F_name_generic_template

    Defaults to '{underscore_name'} - method1

F_name_impl_method_template

    {lower_class}_{underscore_name}{function_suffix}

F_name_impl_function_template

    {underscore_name}{function_suffix}

F_name_method_template

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

   Name of the Fortran argument which is the derived type
   which represents a C++ class.
   It must not be the same as any of the routines arguments.
   Defaults to ``obj``.

F_result

    The name of the Fortran wrapper's result variable.
    It must not be the same as any of the routines arguments.
    It defaults to *rv*  (return value).

F_derived_member

    The name of the member of the Fortran derived type which
    wraps a C++ class.  It will contain a ``type(C_PTR)`` which
    points to the C++ instance.
    Defaults to *voidptr*.


Top Level
---------

copyright

   A list of lines to add to the top of each generate file.

splicers

   A dictionary mapping file suffix to a list of splicer files
   to read.

types

   A dictionary of user define types.
   Each type is a dictionary for members describing how to
   map a type between languages.

patterns:

   Code blocks to insert into generated code.

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

Types
-----

Types describe how to handle arguments from Fortran to C to C++.  Then
how to convert return values from C++ to C to Fortran.

Since Fortran 2003 (ISO/IEC 1539-1:2004(E)) there is a standardized
way to generate procedure and derived-type declarations and global
variables which are interoperable with C (ISO/IEC 9899:1999). The
bind(C) attribute has been added to inform the compiler that a symbol
shall be interoperable with C; also, some constraints are added. Note,
however, that not all C features have a Fortran equivalent or vice
versa. For instance, neither C's unsigned integers nor C's functions
with variable number of arguments have an equivalent in
Fortran. [#f1]_


.. list from util.py class Typedef

base

    Base type.
    For example, string and string_from_buffer both have a 
    base time of *string*.
    Defaults to *unknown*

forward

    Forward declaration.
    Defaults to *None*.

typedef

    Initialize from existing type
    Defaults to *None*.

cpp_type

    Name of type in C++.
    Defaults to *None*.

cpp_to_c

    Expression to convert from C++ to C.
    Defaults to *{var}*.

cpp_header

    Name of C++ header file required for implementation.
    For example, if cpp_to_c was a function.
    Defaults to *None*.

c_type

    name of type in C.
    Defaults to *None*.

c_header

    Name of C header file required for type.
    Defaults to *None*.

c_to_cpp

    Expression to convert from C to C++.
    Defaults to *{var}*.

c_fortran

    Expression to convert from C to Fortran.
    Defaults to *None*.

c_argdecl

    List of argument declarations for C wrapper, *None*=match declaration.
    Used with string_from_buffer .
    Defaults to *None*.

c_intent_in

    Code to add for argument with intent(IN).
    Can be used to convert types or copy-in semantics.
    For example, ``char *`` to ``std::string``.

c_intent_in_trim

    Code to add for argument with intent(IN) and len_trim attribute 
    For example, ''char *, int`` into ``std::string``

c_intent_out

    Code to add after call when ``intent(OUT)`` or ``intent(INOUT)``.
    Used to implement copy-out semantics.

c_pre_call

    Statement to execute before call.

c_post_call

    Statement to execute after call.
    Can be use to cleanup after *c_pre_call*
    or to coerce the return value.
    Defaults to *None*.

c_return_code

    Fortran code used to call function and assign the return value.
    Defaults to *None*.

f_c_args

    List of argument names to F_C routine.
    Defaults to *None*.

f_c_argdecl

    List of declarations to F_C routine.
    By default, only a single argument is passed for each dummy argument.
    Defaults to *None*.

f_type

    Name of type in Fortran.
    Defaults to *None*.

f_derived_type

    Fortran derived type name.
    Defaults to *None* i.e. use C++ class name.

f_args

    Arguments in the Fortran wrapper to pass to the C function.
    This can pass multiple arguments to C for a single
    argument to the wrapper; for example, an address and length
    for a ``character(*)`` argument.
    Or it may be intermediate values.
    For example, a Fortran character variable can be converted
    to a ``NULL`` terminated string with
    ``trim({var}) // C_NULL_CHAR``.
    Defaults to *None*  i.e. pass argument unchanged.

f_argsdecl

    A list of declarations needed by *f_args*, *f_pre_call* or
    *f_post_call*.
    Defaults to *None* i.e. no additional declarations.

f_module

    Fortran modules needed for type  (dictionary).
    Defaults to *None*.

f_return_code

    Fortran code used to call function and assign the return value.
    Defaults to *None*.

.. f_kind

..    Fortran kind of type.
..    Defaults to *None*.

f_cast

    Expression to convert Fortran type to C type.
    This is used when creating a Fortran generic functions which
    accept several type but call a single C function which expects
    a specific type.
    For example, type ``int`` is defined as ``int({var}, C_INT)``.
    This expression converts *var* to a ``integer(C_INT)``.
    Defaults to *{var}*  i.e. no conversion.

f_use_tmp

    If *true*, pass {tmp_var} to C routine instead of {var}.
    This can be used with *f_pre_call* to convert Fortran values
    to values.  For example, to cast or map values.
    Defaults to *False*.

f_pre_call

    Statement to execute before call, often to coerce types
    when *f_cast* cannot be used.  If this involves the temporary
    variable then *f_use_tmp* should be set to *True*.
    Defaults to *None*.

f_post_call

    Statement to execute after call.
    Can be use to cleanup after *f_pre_call*
    or to coerce the return value.
    Defaults to *None*.

..  XXX - maybe later.  For not in wrapping routines
..         f_attr_len_trim = None,
..         f_attr_len = None,
..         f_attr_size = None,

result_as_arg

    Override fields when result should be treated as an argument.
    Defaults to *None*.

PY_format

    'format unit' for PyArg_Parse.
    Defaults to *O*

PY_PyTypeObject

    Variable name of PyTypeObject instance.
    Defaults to *None*.

PY_PyObject

    Typedef name of PyObject instance.
    Defaults to *None*.

PY_ctor

    Expression to create object.
    ex. PyBool_FromLong({rv})
    Defaults to *None*.

PY_to_object

    PyBuild - object = converter(address).
    Defaults to *None*.

PY_from_object

    PyArg_Parse - status = converter(object, address).
    Defaults to *None*.

PY_post_parse

   Used if PY_PyTypeObject is set.
   A format expression to convert a *PyObject* into the type.
   Ex. ``{var} = PyObject_IsTrue({var_obj});``

Format dictionary for Type fields

  * var - name of variable, defaults to argument name.
  * tmp_var - temporary variable.  defaults to *tmp_{var}*.
  * result_arg - name of result variable from *F_string_result_as_arg*.
  * F_result - name of result variable
  * F_C_name - name of BIND(C) interface
  * F_arg_c_call
  * F_arg_c_call_tab
  * F_arguments


arg_f_decl._f_decl(arg)

Example for each type::

   subroutine name({var})
       {f_argsdecl}

       ! arguments
       foreach argument:
          F_arg_c_call += f_args or f_cast or '{var}'

       {f_pre_call}
       {f_return_code}     ! call C code
       {f_post_call}



Predefined types

  * void
  * int
  * long
  * size_t
  * float
  * double
  * bool
  * string
  * string_from_buffer


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

default_arg_suffix

   A list of suffixes to apply to C and Fortran functions generated when
   wrapping a C++ function with default arguments.  The first entry is for
   the function with the fewest arguments and the final entry should be for
   all of the arguments.

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

F_name_method

    The name of the *F_name_impl* subprogram when used as a
    type procedure.
    Defaults to option *F_name_method_template*.

F_name_generic

..    method1
    Defaults to option *F_name_generic_template*.



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

len

   An expression for the length of string result variable.
   If not set then the function will be called to compute the string
   result and len will be computed using ``strlen``.
   The function is then called again to fill in the result variable.
 
len_trim

   For a string argument, pass the string address and the result of
   len_trim.

Doxygen
-------

Used to insert directives for doxygen.

brief

   Brief description.

description

   Full description.

return

   Description of return value.


Splicers
--------

Describe splicers.



.. rubric:: Footnotes

.. [#f1] https://gcc.gnu.org/onlinedocs/gfortran/Interoperability-with-C.html

