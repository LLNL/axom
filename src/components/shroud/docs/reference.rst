Reference
=========

Command Line Options
--------------------

help
       Show this help message and exit

version
       Show program's version number and exit

outdir OUTDIR
       Directory for output files

outdir-c-fortran OUTDIR_C_FORTRAN
       Directory for C/Fortran wrapper output files, overrides --outdir

outdir-python OUTDIR_PYTHON
       Directory for Python wrapper output files, overrides --outdir

logdir LOGDIR
       Directory for log files

cfiles CFILES
       Output file with list of C and C++ files created

ffiles FFILES
       Output file with list of Fortran created

path PATH
       Colon delimited paths to search for splicer files, may
       be supplied multiple times to create path



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

F_C_prefix
  Prefix added to name of generated Fortran interface for C routines.
  Defaults to **c_**.

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

wrap_lua
  If *true*, create Lua wrappers.
  Defaults to *false*.


Option Templates
^^^^^^^^^^^^^^^^

Templates are set in options then expanded to assign to the format 
dictionary.

C_name_template
    ``{C_prefix}{class_name}_{underscore_name}{function_suffix}``

F_C_name_template
    ``{F_C_prefix}{class_name}{underscore_name}{function_suffix}``

F_name_generic_template
    ``{underscore_name}``

F_name_impl_template
    ``{class_name}{underscore_name}{function_suffix}``

F_name_method_template
    ``{underscore_name}{function_suffix}``

PY_name_impl
    PY_class1_method1




C_header_filename_library_template
   ``wrap{library}.h``

C_impl_filename_library_template
    ``wrap{library}.cpp``

C_header_filename_class_template
    ``wrap{cpp_class}.h``

C_impl_filename_class_template
    ``wrap{cpp_class}.cpp``


F_module_name_library_template
    ``{lower_library}_mod``

F_impl_filename_library_template
    ``wrapf{lower_library}.f``

F_module_name_class_template
    ``{class_lower}_mod``

F_impl_filename_class_template
    ``wrapf{cpp_class}.f``

F_name_impl_template
    ``{name_class}{underscore_name}{function_suffix}``


LUA_module_filename_template
    ``lua{library}module.cpp``

LUA_header_filename_template
    ``lua{library}module.hpp``

LUA_userdata_type_template
    ``{LUA_prefix}{cpp_class}_Type``

LUA_userdata_member_template
    Name of pointer to class instance in userdata.
    ``self``

LUA_class_reg_template
    Name of `luaL_Reg` array of function names.
    ``{LUA_prefix}{cpp_class}_Reg``

LUA_metadata_template
    Name of metatable for a class.
    ``{cpp_class}.metatable``

LUA_ctor_name_template
    Name of constructor for a class.
    Added to the library's table.
    ``{cpp_class}``

LUA_name_impl_template
    ``{LUA_prefix}{class_name}{underscore_name}{function_suffix}``





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


Top Level Fields
----------------

copyright
   A list of lines to add to the top of each generate file.

splicers
   A dictionary mapping file suffix to a list of splicer files
   to read.

types
   A dictionary of user define types.
   Each type is a dictionary for members describing how to
   map a type between languages.

patterns
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

Types Dictionary
----------------

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
    Defaults to *{cpp_var}*.

cpp_header
    Name of C++ header file required for implementation.
    For example, if cpp_to_c was a function.
    Defaults to *None*.

c_type
    name of type in C.
    Defaults to *None*.

c_header
    Name of C header file required for type.
    This file is included in the interface header.
    Defaults to *None*.

c_to_cpp
    Expression to convert from C to C++.
    Defaults to *{c_var}*.

c_fortran
    Expression to convert from C to Fortran.
    Defaults to *None*.

c_statements
    A nested dictionary of code template to add.
    The first layer is *intent_in*, *intent_out*, and *result*.
    The second layer is *pre_call*, *pre_call_trim*, *post_call*, *cpp_header*.
    The entries are a list of format strings.

    intent_in
        Code to add for argument with intent(IN).
        Can be used to convert types or copy-in semantics.
        For example, ``char *`` to ``std::string``.

    intent_in_trim
        Code to add for argument with intent(IN) and len_trim attribute 
        For example, ``char *, int`` into ``std::string``

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

f_statement
    A nested dictionary of code template to add.
    The first layer is *intent_in*, *intent_out*, and *result*.
    The second layer is *declare*, *pre_call*, and *post_call*
    The entries are a list of format strings.

    declare
        A list of declarations needed by *pre_call* or *f_post_call*.

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
   to chain calls in C++.  This does not translate to C and Fortran.
   Instead make the return type ``void``.



C_name
    Name of the C wrapper function.
    Defaults to option *C_name_template*.

F_C_name
    Name of the Fortran ``BIND(C)`` interface for a C function.
    Defaults to the lower case version of *F_C_name_template*.

..    tut_class1_method1

F_name_impl
    Name of the Fortran implementation function.
    Defaults to option *F_name_impl_template* .

..    class1_method1

F_name_method
    The name of the *F_name_impl* subprogram when used as a
    type procedure.
    Defaults to option *F_name_method_template*.

F_name_generic
    Defaults to option *F_name_generic_template*.

F_name_instance_get
    Name of method to get ``type(C_PTR)`` instance pointer from wrapped class.
    Defaults to *get_instance*.
    If the name is blank, no function is generated.

F_name_instance_set
    Name of method to set ``type(C_PTR)`` instance pointer in wrapped class.
    Defaults to *set_instance*.
    If the name is blank, no function is generated.

Annotations
-----------

a.k.a. attributes

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
   An expression for the length of string result variable.
   If not set then the function will be called to compute the string
   result and len will be computed using ``strlen``.
   The function is then called again to fill in the result variable.
 
len_trim
   For a string argument, pass the string address and the result of
   len_trim.

Doxygen
-------

Used to insert directives for doxygen for a function.

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

