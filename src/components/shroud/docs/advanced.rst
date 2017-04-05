Advanced Usage
==============

Customizing Behavior
--------------------

Fields
^^^^^^

Fields only apply to the type, class or function to which it belongs.
They are not inherited.
For example, *C_name* is a field which is used to explicitly name
a single C wrapper function.  While *C_name_template* is an option which
controls the default value of *C_name*.

Annotations
^^^^^^^^^^^

Annotations or attributes apply to specific arguments or results.
They describe semantic behavior for an argument.

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

What files are created
----------------------

Shroud will create multiple output file which must be compiled with C++ or Fortran compilers.

One C++ file will be created for the library and one file for each C++ class.

By default, Fortran will create one file per class similar to the way C is handled.

If one class makes use of another class, it is necessary to put all of the class
into a single file using the *F_module_per_class* option.
Each Fortran file will only contain one module to make it easier to create makefile
dependencies using pattern rules::

    %.o %.mod : %.f


How Names are Computed
----------------------

Shroud attempts to provide user control of names while providing reasonable defaults.
Names are controlled by a format string or can be specified explicitly.

Format strings contain “replacement fields” surrounded by curly braces
{}. Anything that is not contained in braces is considered literal
text, which is copied unchanged to the output. If you need to include
a brace character in the literal text, it can be escaped by doubling:
{{ and }}. [Python_Format]_

The replacement fields are defined by the format dictionary.  Shroud
defines values which may be used.

Library name - Updated after reading YAML file.
   * library - The value of **field** *library*.
   * library_lower - Lowercase version of *library*.
   * library_upper - Uppercase version of *library*.

Class name - Updated before processing each class.
   * cpp_class - The name of the C++ class from the YAML input file.
   * class_lower - Lowercase version of *cpp_class*.
   * class_upper - Uppercase version of *cpp_class*.
   * class_name  - Variable which may be used in creating function names.
                   Defaults to evaluation of *class_name_template*.
                   Outside of a class, set to empty string.
   * C_prefix - Prefix for C wrapper functions.
     Defaults to first three letters of *library*.
     Set from **options**.
   * F_C_prefix - Prefix for Fortran name for C wrapper.  Defaults to ``c_``.
     Set from **options**.

Function name - Updated before processing each function or method.
   * function_name - Name of function in the YAML file.
   * underscore_name - *function_name* converted from CamelCase to snake_case.
   * function_suffix - Suffix append to name.  Used to differentiate overloaded functions.
     Defaults to a sequence number (e.g. `_0`, `_1`, ...) but can be set
     by using the function field *function_suffix*.
     Mulitple suffixes may be applied.



+------------------------+---------------------------------+------------------+
| Description            | Option                          | Override         |
+========================+=================================+==================+
| C wrapper              | *C_name_template*               | *C_name*         |
| implementation         |                                 |                  |
+------------------------+---------------------------------+------------------+
| Fortran BIND(C)        | *F_C_name_template*             | *F_C_name*       |
| interface              |                                 |                  |
+------------------------+---------------------------------+------------------+
| Fortran wrapper        | *F_name_impl_template*          | *F_name_impl*    |
| implementation         |                                 |                  |
+------------------------+---------------------------------+------------------+
| Fortran method         | *F_name_method_template*        | *F_name_method*  |
+------------------------+---------------------------------+------------------+
| Fortran generic name   | *F_name_generic_template*       | *F_name_generic* |
+------------------------+---------------------------------+------------------+

example
^^^^^^^

For example, a library has an ``initialize`` function which is
in a namespace.  In C++ it is called as::

  #include "library.hpp"

  library::initialize()

By default this will be a function in a Fortran module and 
can be called as::

  use library

  call initialize

Since ``initialize`` is a rather common name for a function, it may 
be desirable to rename the Fortran wrapper to something more specific.
The name of the Fortran implementation wrapper can be changed
by setting *F_name_impl*::

  options:
    library: library
    namespace: library

  function:
  -  decl: void initialize
     F_name_impl: library_initialize

To rename all functions, set the template in the toplevel *options*::     

  options:
    library: library
    namespace: library
    F_name_impl_function_template:
      "{library}_{underscore_name}{function_suffix}"

  function:
  -  decl: void initialize

The alternative is to have the user rename it at the point
of use if there is a conflict with other modules::

   use library, library_initialize => initialize

   call library_initialize


How Functions are Generated
---------------------------

This section show the format templates which are used to create code.
The names in curly parens are user settable values. 

The C wrapper code::

    struct s_{C_type_name};
    typedef struct s_{C_type_name} {C_type_name};

C implementation::

    {C_return_type} {C_name}({C_prototype})
    {
        {C_const}{cpp_class} *{C_this}obj =
            static_cast<{C_const}{cpp_class} *>(
                static_cast<{C_const}void *>({C_this}));

        {rv_decl} = {CPP_this_call}{method_name}{CPP_template}({C_call_list});
        // pre_call
        {C_code}
        // post-call
        // return_line
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
        ! pre-call
        {F_code}
      end subroutine {F_name_impl}

    end module {F_module_name}


Code
^^^^

Each argument may add code to the *pre_call* and *post_call* sections.

Return value.

Example code::

    {C_return_type} {C_name}({C_prototype})
    {
        {C_const}{cpp_class} *{C_this}obj = new {cpp_class}({C_call_list});
        {C_code}
        return static_cast<AA_exclass1 *>(static_cast<void *>(selfobj));
    }



        ExClass1 *selfobj = new ExClass1(name);


Annotations may change how the code is generated.
The *constructor* attribute will use the `new` C++ keyword and
*destructor* will use `delete` in the *C_code*.


Header Files
^^^^^^^^^^^^

The header files for the library are included by the generated C++ source files.

The library source file will include the global *cpp_header* field.
Each class source file will include the class *cpp_header* field unless it is blank.
In that case the global *cpp_header* field will be used.

To include a file in the implementation list it in the global or class options::

    cpp_header: global_header.hpp

    classes:
    -  name: Class1
       cpp_header: class_header.hpp

    types:
       CustomType:
          typedef: int
          c_header:  type_header.h
          cpp_header : type_header.hpp


The *c_header* field will be added to the header file of contains functions
which reference the type.
This is used for files which are not part of the library but which contain code
which helps map C++ constants to C constants

.. FILL IN MORE

Namespace
---------

Each library or class can be associated with a namespace::

    namespace one {
    namespace two {
       void function();

       namespace three {
         class Class1 {
         };
       }

       class Class2 {
       };
    }
    }

The YAML file would look like::

    namespace: one two

    classes:
    -  Class1
       cpp_header: one two three
    -  Class2


Local Variable
^^^^^^^^^^^^^^

*SH_* prefix on local variables.

Results are named from *fmt.rv*.

Fortran option F_result.


Character Type
--------------

Fortran, C, and C++ all have their own semantics for character variables.

  * Fortran ``character`` variables know their length and are blank filled
  * C ``char *`` variables are assumed to be ``NULL`` terminated.
  * C++ ``std::string`` know their own length and are ``NULL`` terminated.

It is not sufficient to pass an address between Fortran and C++ like
it is with other native types.  In order to get ideomatic behavior in
the Fortran wrappers it is often necessary to copy the values.  This
is to account for blank filled vs ``NULL`` terminated.  It also helps
support ``const`` vs non-``const`` strings.

A C 'bufferify' wrapper is created which accepts the address of the
Fortran character variable with a ``int`` argument for the declared
length of the variable (``len``) and/or a ``int`` argument for the
length with blanks trimmed off (``len_trim``).
The wrapper then uses these arguments to create a ``NULL`` terminated string
or a std::string instance.

Character Arguments
^^^^^^^^^^^^^^^^^^^

When an argument has intent *out*, then *len* attribute is added.
This allows the wrapper routine to know how much space as available for the output string.

When the argument has intent *in*, then the *len_trim* attribute is added to the *bufferify*
wrapper only.  The non-bufferify version will use ``strlen`` to compute the length of data.

Character Function
^^^^^^^^^^^^^^^^^^

.. This stuff was moved here from the tutorial and should be cleaned up

This attribute marks the routine as Fortran ``pure`` meaning there are
no side effects.  This is necessary because the function will be
called twice.  Once to compute the length of the result and once to
return the result.

The length of result variable ``rv`` is computed by calling the
function.  Once the result is declared, ``tut_function4a`` is called
which returns a ``type(C_PTR)``.  This result is dereferenced by
``fstr`` and copied into ``rv``.


.. XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

It is possible to avoid calling the C++ function twice by passing in
another argument to hold the result.  It would be up to the caller to
ensure it is long enough.  This is done by setting the option
**F_string_result_as_arg** to true.  Like all options, it may also be
set in the global **options** and it will apply to all functions::

.. update code examples from current output






    - decl: const std::string& Function4b(
        const std::string& arg1,
        const std::string& arg2)
      options:
        F_string_result_as_arg: output

The generated Fortran wrapper::

    subroutine function4b(arg1, arg2, output)
        use iso_c_binding, only : C_INT
        implicit none
        character(*), intent(IN) :: arg1
        character(*), intent(IN) :: arg2
        character(*), intent(OUT) :: output
        rv = c_function4b_bufferify(  &
            arg1,  &
            len_trim(arg1),  &
            arg2,  &
            len_trim(arg2),
            output,  &
            len(output))
    end subroutine function4b

The generated C wrapper::

    void TUT_function4b_bufferify(const char * arg1, int Larg1,
                                  const char * arg2, int Larg2,
                                  char * output, int Loutput) {
        const std::string SH_arg1(arg1, Larg1);
        const std::string SH_arg2(arg2, Larg2);
        const std::string & rv = Function4b(SH_arg1, SH_arg2);
        asctoolkit::shroud::FccCopy(output, Loutput, rv.c_str());
        return;
    }


 ``FccCopy`` will copy the result into ``output`` and blank fill.


.. char **


Complex Type
------------


Derived Types
-------------



* chained function calls


splicers
--------


.. [Python_Format] https://docs.python.org/2/library/string.html#format-string-syntax




