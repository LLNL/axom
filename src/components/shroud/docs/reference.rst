Reference
=========

Command Line Options
--------------------

help
       Show this help message and exit.

version
       Show program's version number and exit.

outdir OUTDIR
       Directory for output files.
       Defaults to current directory.

outdir-c-fortran OUTDIR_C_FORTRAN
       Directory for C/Fortran wrapper output files, overrides *--outdir*.

outdir-python OUTDIR_PYTHON
       Directory for Python wrapper output files, overrides *--outdir*.

outdir-lua OUTDIR_LUA
       Directory for Lua wrapper output files, overrides *--outdir*.

outdir-yaml OUTDIR_YAML
       Directory for YAML output files, overrides *--outdir*.

logdir LOGDIR
       Directory for log files.
       Defaults to current directory.

cfiles CFILES
       Output file with list of C and C++ files created.

ffiles FFILES
       Output file with list of Fortran created.

path PATH
       Colon delimited paths to search for splicer files, may
       be supplied multiple times to append to path.

sitedir
       Return the installation directory of shroud and exit.
       This path can be used to find cmake/SetupShroud.cmake.


Format Dictionary
-----------------

Each scope has its own format dictionary.  If a value is not found in
the dictionary, then the parent dictionaries will be recursively
searched.

Library
^^^^^^^

C_header_filename
    Name of generated header file for the library.
    Defaulted from expansion of option *C_header_filename_library_template*.

C_impl_file
    Name of generated C++ implementation file for the library.
    Defaulted from expansion of option *C_impl_filename_library_template*.

C_result
    TODO

F_module_name
    Name of module for Fortran interface for the library.
    Defaulted from expansion of option *F_module_name_library_template*.

F_impl_filename
    Name of generated Fortran implemention file for the library.
    Defaulted from expansion of option *F_impl_filename_library_template*.
    If option *F_module_per_class* is false, then all derived types
    generated for each class will also be in this file.

F_result
    TODO

library
    The value of global **field** *library*.

library_lower
    Lowercase version of *library*.

library_upper
    Uppercase version of *library*.

namespace_scope
    The values in field **namespace** delimited with ``::``.


Class
^^^^^

C_header_filename
    Name of generated header file for the class.
    Defaulted from expansion of option *C_header_filename_class_template*.

C_impl_file
    Name of generated C++ implementation file for the library.
    Defaulted from expansion of option *C_impl_filename_class_template*.

F_module_name
    Name of module for Fortran interface for the library.
    Defaulted from expansion of option *F_module_name_class_template*.
    Only defined if *F_module_per_class* is true.

F_impl_filename
    Name of generated Fortran implemention file for the library.
    Defaulted from expansion of option *F_impl_filename_class_template*.
    Only defined if *F_module_per_class* is true.

cpp_class
    The name of the C++ class from the YAML input file.

class_lower
    Lowercase version of *cpp_class*.

class_upper
    Uppercase version of *cpp_class*.

class_prefix
    Variable which may be used in creating function names.
    Defaults to evaluation of *class_prefix_template*.
    Outside of a class, set to empty string.

C_prefix
    Prefix for C wrapper functions.
    Set from **options**.
    If no option is set then the first three letters
    of *library_upper* are used.

F_C_prefix
    Prefix for Fortran name for C wrapper.  Defaults to ``c_``.
    Set from **options** and defaults to `c_`



Function
^^^^^^^^

C_call_list
    Comma delimited list of function arguments.

C_post_call
    Statements added after the call to the function.
    Used to convert result and/or ``intent(OUT)`` arguments to C types.

C_pre_call
    Statements added before the call to the function.
    Used to convert C types to C++ types.

C_prototype
    C prototype for the function.
    This will include any arguments required by annotations or options,
    such as length or **F_string_result_as_arg**.  

C_return_type
    Return type of the function.
    If the **return_this** field is true, then *C_return_type* is set to ``void``.

CPP_template
    The template component of the function declaration.
    ``<{type}>``

CPP_this_call
    How to call the function.
    ``{CPP_this}->`` for instance methods and blank for library functions.

F_arg_c_call
    Comma delimited arguments to call C function from Fortran.

F_arg_c_call_tab
    Tab delimited version *F_arg_c_call*.
    Used to avoid long lines.

F_arguments
    Set from option *F_arguments* or generated from YAML decl.

F_C_arguments
    Argument names to the ``bind(C)`` interface for the subprogram.

F_C_call
    The name of the C function to call.  Usually *F_C_name*, but it may
    be different if calling a generated routine.
    This can be done for functions with string arguments.

F_C_name
    The name of the ``bind(C)`` interface function.

F_C_pure_clause
    TODO

F_C_result_clause
    Result clause for the ``bind(C)`` interface.

F_C_subprogram
    ``subroutine`` or ``function``.

F_pure_clause
    For non-void function, ``pure`` if the *pure* annotation is added or 
    the function is ``const`` and all arguments are ``intent(in)``.

F_name_method
    Evaluation of *F_name_method_template*.

F_name_impl
    Evaluate of *F_name_impl_template*.

F_result_clause
    `` result({F_result})`` for functions.
    Blank for subroutines.

function_name
    Name of function in the YAML file.

underscore_name
    *function_name* converted from CamelCase to snake_case.

function_suffix
    Suffix append to name.  Used to differentiate overloaded functions.
    Defaults to a sequence number (e.g. `_0`, `_1`, ...) but can be set
    by using the function field *function_suffix*.
    Mulitple suffixes may be applied.

C_rv_decl
    Declaration of return value for function.

Argument
^^^^^^^^

C_const
    ``const`` if argument has the *const* attribute.

c_var
    The C name of the argument.

c_var_len
    Function argument generated from the *len* annotation.
    Set from option **C_var_len_template**.

c_var_trim
    Function argument generated from the *len_trim* annotation.
    Set from option **C_var_trim_template**.

cpp_type
    The C++ type of the argument.

cpp_var
    Name of the C++ variable.

cpp_val
    Evaluation of cpp_to_c for the arguments typedef.

f_var
    Fortran variable name for argument.

c_ptr
    `` * `` if argument is a pointer.

len_var
    TODO

Global Fields
-------------

C_header_filename
   Output file name for header for  wrapper routines.
   Defaults to expansion of option *C_header_filename_library_template*.

C_impl_filename
   Output file name for implementation of wrapper routines.
   Defaults to expansion of option *C_impl_filename_library_template*.

copyright
   A list of lines to add to the top of each generate file.
   Do not include any language specific comment characters since
   Shroud will add the appropriate comment delimiters for each language.

cpp_header
  C++ header file name which will be included in the implementation file.

F_module_name
   Name of Fortran module for this class.
   Defaults to option *F_module_name_library_template*.

F_impl_filename
   Name of Fortran file for functions.
   Defaults to option *F_impl_name_library_template*.

library
  The name of the library.
  Used to name output files and modules.
  The first three letters are used as the default for **C_prefix** option.
  Defaults to *default_library*.
  Each YAML file is intended to wrap a single library.

namespace
  Blank delimited list of namespaces for **cpp_header**.
  The namespaces will be nested.

patterns
   Code blocks to insert into generated code.

splicers
   A dictionary mapping file suffix to a list of splicer files
   to read.

types
   A dictionary of user define types.
   Each type is a dictionary for members describing how to
   map a type between languages.

Options
-------

debug
  Print additional comments in generated files that may 
  be useful for debugging.
  Defaults to *false*.

C_bufferify_suffix
  Suffix appended to generated routine which pass strings as buffers
  with explicit lengths.
  Defaults to *_bufferify*

C_prefix
  Prefix added to name of generated C routines.
  The prefix helps to ensure unique global names.

C_proto_type
   XXX  override prototype of generated C function

C_result
    The name of the C wrapper's result variable.
    It must not be the same as any of the routines arguments.
    It defaults to *SH_rv*  (Shroud return value).

C_return_type
   XXX   override return type of function

C_string_result_as_arg
  The name of the output argument for string results.
  Function which return ``char`` or ``std::string`` values return
  the result in an additional argument in the C wrapper.

C_this
    Name of the C object argument.  Defauls to ``self``.
    It may be necessary to set this if it conflicts with an argument name.

C_var_len_template
    Format for variable created with *len* annotation.
    Default ``N{c_var}``

C_var_trim_template
    Format for variable created with *len_trim* annotation.
    Default ``L{c_var}``

CPP_this
    Name of the C++ object pointer set from the *C_this* argument.
    Defauls to ``SH_this``.


F_C_prefix
  Prefix added to name of generated Fortran interface for C routines.
  Defaults to **c_**.

F_derived_member
    The name of the member of the Fortran derived type which
    wraps a C++ class.  It will contain a ``type(C_PTR)`` which
    points to the C++ instance.
    Defaults to *voidptr*.

F_this
   Name of the Fortran argument which is the derived type
   which represents a C++ class.
   It must not be the same as any of the routines arguments.
   Defaults to ``obj``.

F_result
    The name of the Fortran wrapper's result variable.
    It must not be the same as any of the routines arguments.
    It defaults to *SH_rv*  (Shroud return value).

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


LUA_result
    The name of the Lua wrapper's result variable.
    It defaults to *rv*  (return value).

PY_result
    The name of the Python wrapper's result variable.
    It defaults to *rv*  (return value).

show_splicer_comments
    If ``true`` show comments which delineate the splicer blocks;
    else, do not show the comments.
    Only the global level option is used.

wrap_c
  If *true*, create C wrappers.
  Defaults to *true*.

wrap_fortran
  If *true*, create Fortran wrappers.
  Defaults to *true*.

wrap_python
  If *true*, create Python wrappers.
  Defaults to *false*.

wrap_lua
  If *true*, create Lua wrappers.
  Defaults to *false*.


Option Templates
^^^^^^^^^^^^^^^^

Templates are set in options then expanded to assign to the format 
dictionary.

C_header_filename_class_template
    ``wrap{cpp_class}.h``

C_header_filename_library_template
   ``wrap{library}.h``

C_impl_filename_class_template
    ``wrap{cpp_class}.cpp``

C_impl_filename_library_template
    ``wrap{library}.cpp``

C_name_template
    ``{C_prefix}{class_prefix}{underscore_name}{function_suffix}``

class_prefix_template
    Class component for function names.
    Will be blank if the function is not in a class.
    ``{class_lower}_``

F_C_name_template
    ``{F_C_prefix}{class_prefix}{underscore_name}{function_suffix}``

F_name_generic_template
    ``{underscore_name}``

F_impl_filename_class_template
    ``wrapf{cpp_class}.f``

F_impl_filename_library_template
    ``wrapf{library_lower}.f``

F_name_impl_template
    ``{name_class}{underscore_name}{function_suffix}``

F_module_name_class_template
    ``{class_lower}_mod``

F_module_name_library_template
    ``{library_lower}_mod``

F_name_impl_template
    ``{class_prefix}{underscore_name}{function_suffix}``

F_name_function_template
    ``{underscore_name}{function_suffix}``


LUA_class_reg_template
    Name of `luaL_Reg` array of function names for a class.
    ``{LUA_prefix}{cpp_class}_Reg``

LUA_ctor_name_template
    Name of constructor for a class.
    Added to the library's table.
    ``{cpp_class}``

LUA_header_filename_template
    ``lua{library}module.hpp``

LUA_metadata_template
    Name of metatable for a class.
    ``{cpp_class}.metatable``

LUA_module_filename_template
    ``lua{library}module.cpp``

LUA_module_name
    Name of Lua module for library.
    ``{library_lower}``

LUA_module_reg_template
    Name of `luaL_Reg` array of function names for a library.
    ``{LUA_prefix}{library}_Reg``

LUA_name_impl_template
    Name of implementation function.
    All overloaded function use the same Lua wrapper so 
    *function_suffix* is not needed.
    ``{LUA_prefix}{class_prefix}{underscore_name}``

LUA_name_template
    Name of function as know by Lua.
    All overloaded function use the same Lua wrapper so 
    *function_suffix* is not needed.
    ``{function_name}``

LUA_userdata_type_template
    ``{LUA_prefix}{cpp_class}_Type``

LUA_userdata_member_template
    Name of pointer to class instance in userdata.
    ``self``


PY_name_impl
    PY_class1_method1


Types Map
---------

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
    Defaults to *{cpp_var}*.  i.e. no conversion required.

cpp_header
    Name of C++ header file required for implementation.
    For example, if cpp_to_c was a function.
    Defaults to *None*.

cpp_local_var
    If true then a local variable will be created instead of passing the argument
    directly to the function.
    The variable will be assigned a value using *c_to_cpp*.
    If *c_to_cpp* is a large expression it is sometimes convient to have a local variable
    for debugging purposes.
    It can also be used to create cleaner code when *c_to_cpp* will generate a very long statement.
    When *c_to_cpp* is not sufficient to assign a value, *c_statements* can be used to 
    add multiple statements into the wrapper.  *c_statements* and *cpp_local_var* cannot
    be used together.

..  {C_const}{cpp_type}{ptr} = c_to_cpp ;

c_type
    name of type in C.
    Defaults to *None*.

c_header
    Name of C header file required for type.
    This file is included in the interface header.
    Defaults to *None*.

c_to_cpp
    Expression to convert from C to C++.
    Defaults to *{c_var}*.  i.e. no conversion required.

c_statements
    A nested dictionary of code template to add.
    The first layer is *intent_in*, *intent_out*, *result*,
    *intent_in_buf*, *intent_out_buf*, and *result_buf*.
    The second layer is *pre_call*, *pre_call_buf*, *post_call*, *cpp_header*.
    The entries are a list of format strings.

    intent_in
        Code to add for argument with intent(IN).
        Can be used to convert types or copy-in semantics.
        For example, ``char *`` to ``std::string``.

    intent_out
        Code to add after call when ``intent(OUT)`` or ``intent(INOUT)``.
        Used to implement copy-out semantics.

    result
        Code to use when passing result as an argument.

        cpp_header
           string of blank delimited header names

        cpp_local_var
           True if a local C++ variable is created.
           This is the case when C and C++ are not directly compatible.
           Usually a C++ constructor is involved.
           This sets *cpp_var* is set to ``SH_{c_var}``.

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

f_c_module
    Fortran modules needed for type in the interface.
    A dictionary keyed on the module name with the value being a list of symbols.
    Similar to **f_module**.
    Defaults to *None*.

f_c_type
    Type declaration for ``bind(C)`` interface.
    Defaults to *None* which will then use *f_type*.

f_type
    Name of type in Fortran.
    Defaults to *None*.

f_derived_type
    Fortran derived type name.
    Defaults to *None* which will use the C++ class name
    for the Fortran derived type name.

.. f_args
    Arguments in the Fortran wrapper to pass to the C function.
    This can pass multiple arguments to C for a single
    argument to the wrapper; for example, an address and length
    for a ``character(*)`` argument.
    Or it may be intermediate values.
    For example, a Fortran character variable can be converted
    to a ``NULL`` terminated string with
    ``trim({var}) // C_NULL_CHAR``.
    Defaults to *None*  i.e. pass argument unchanged.

f_module
    Fortran modules needed for type in the implementation wrapper.
    A dictionary keyed on the module name with the value being a list of symbols.
    Defaults to *None*.::

        f_module:
           iso_c_binding:
             - C_INT

f_return_code
    Fortran code used to call function and assign the return value.
    Defaults to *None*.

f_cast
    Expression to convert Fortran type to C type.
    This is used when creating a Fortran generic functions which
    accept several type but call a single C function which expects
    a specific type.
    For example, type ``int`` is defined as ``int({f_var}, C_INT)``.
    This expression converts *f_var* to a ``integer(C_INT)``.
    Defaults to *{f_var}*  i.e. no conversion.

..  See tutorial function9 for example.  f_cast is only used if the types are different.

f_to_c
    Expression to convert Fortran type to C type.
    If this field is set, it will be used before f_cast.
    Defaults to *None*.

f_statement
    A nested dictionary of code template to add.
    The first layer is *intent_in*, *intent_out*, and *result*.
    The second layer is *declare*, *pre_call*, and *post_call*
    The entries are a list of format strings.

    c_local_var
        If true, generate a local variable using the C declaration for the argument.
        This variable can be used by the pre_call and post_call statements.
        A single declaration will be added even if with ``intent(inout)``.

    declare
        A list of declarations needed by *pre_call* or *f_post_call*.
        Usually a *c_local_var* is sufficient.
        If both *pre_call* and *post_call* are specified then both *declare*
        clause will be added and thus should not declare the same variable.

    pre_call
        Statement to execute before call, often to coerce types
        when *f_cast* cannot be used.

    post_call
        Statement to execute after call.
        Can be use to cleanup after *f_pre_call*
        or to coerce the return value.

    need_wrapper
        If true, the fortran wrapper will always be created.
        This is useful then a function assignment is needed to do a type coercision.

..  XXX - maybe later.  For not in wrapping routines
..         f_attr_len_trim = None,
..         f_attr_len = None,
..         f_attr_size = None,

f_helper
    Additional code to add into the module for helper functions.

    private
       List of names which should be PRIVATE to the module

    interface
       Code to add to the non-executable part of the module.

    source
       Code to add in the CONTAINS section of the module.

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

py_statement
    A nested dictionary of code template to add.
    The first layer is *intent_in*, *intent_out*, and *result*.
    The entries are a list of format strings.

..    declare
        A list of declarations needed by *pre_call* or *f_post_call*.

    post_parse
        Statements to execute after the call to ``PyArg_ParseTupleAndKeywords``.
        Used to convert C values into C++ values.
	Ex. ``{var} = PyObject_IsTrue({var_obj});``

    ctor
        Statements to create a Python object.
	Must ensure that ``py_var = cpp_var`` in some form.

..    post_call
        Statement to execute after call.
        Can be use to cleanup after *f_pre_call*
        or to coerce the return value.

        cpp_local_var
           True if a local C++ variable is created.
           This is the case when C and C++ are not directly compatible.
           Usually a C++ constructor is involved.




Format dictionary for Type fields
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


Class Fields
------------

cpp_header
  C++ header file name which will be included in the implementation file.
  If unset then the global *cpp_header* will be used.

C_header_filename
   Output file name for header for  wrapper routines.
   Defaults to evaluation of option *C_header_filename_class_template*.

C_impl_filename
   Output file name for implementation of wrapper routines.
   Defaults to evaluation of option *C_impl_filename_class_template*.

F_derived_name
   Name of Fortran derived type for this class.
   Defaults to the C++ class name.

F_module_name
   Name of Fortran module for this class.
   Defaults to evaluation of option *F_module_name_class_template*.
   Only used if option *F_module_per_class* is True.

F_impl_filename
   Name of Fortran file for this class.
   Defaults to evaluation of option *F_impl_name_class_template*.
   Only used if option *F_module_per_class* is True.

namespace
  Blank delimited list of namespaces for **cpp_header**.
  The namespaces will be nested.
  If not defined then the global *namespace* will be used.
  If it starts with a ``-`` then no namespace will be used.


Function Fields
---------------

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
   to chain calls in C++.  This ideom does not translate to C and Fortran.
   Instead the *C_return_type* format is set to ``void``.


C_code
    C++ code to use within the splicer block for this function.

C_name
    Name of the C wrapper function.
    Defaults to evaluation of option *C_name_template*.

F_C_name
    Name of the Fortran ``BIND(C)`` interface for a C function.
    Defaults to the lower case version of *F_C_name_template*.

F_code
    Fortran code to use within the splicer block for this function.

..    tut_class1_method1

F_name_impl
    Name of the Fortran implementation function.
    Defaults to evaluation of option *F_name_impl_template* .

..    class1_method1

F_name_function
    The name of the *F_name_impl* subprogram when used as a
    type procedure.
    Defaults to evaluation of option *F_name_function_template*.

F_name_generic
    Defaults to evaluation of option *F_name_generic_template*.

F_name_instance_get
    Name of method to get ``type(C_PTR)`` instance pointer from wrapped class.
    Defaults to *get_instance*.
    If the name is blank, no function is generated.

F_name_instance_set
    Name of method to set ``type(C_PTR)`` instance pointer in wrapped class.
    Defaults to *set_instance*.
    If the name is blank, no function is generated.

LUA_name
    Name of function as known by LUA.
    Defaults to evaluation of option *LUA_name_template*.


Annotations
-----------

Ann annotation can be used to provide semantic information for a function or argument.


.. a.k.a. attributes

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
   Default value for C++ function argument.

len
   For a string argument, pass an additional argument to the
   C wrapper with the result of the Fortran intrinsic ``len``.
   If a value for the attribute is provided it will be the name
   of the extra argument.  If no value is provided then the
   argument name defaults to option *C_var_len_template*.

   When used with a function, it will be the length of the return
   value of the function using the declaration::

     character(kind=C_CHAR, len={c_var_len}) :: {F_result}

len_trim
   For a string argument, pass an additional argument to the
   C wrapper with the result of the Fortran intrinsic ``len_trim``.
   If a value for the attribute is provided it will be the name
   of the extra argument.  If no value is provided then the
   argument name defaults to option *C_var_trim_template*.


Doxygen
-------

Used to insert directives for doxygen for a function.

brief
   Brief description.

description
   Full description.

return
   Description of return value.



.. rubric:: Footnotes

.. [#f1] https://gcc.gnu.org/onlinedocs/gfortran/Interoperability-with-C.html

