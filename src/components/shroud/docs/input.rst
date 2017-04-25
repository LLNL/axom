YAML Input
==========

How code is generated
---------------------

In some sense, Shroud can be thought of as a fancy macro processor.
It takes the function declarations from the YAML file and break them
down into a series of contexts (library, class, function) and defines
a dictionary of format macros of the form key=value.  There are then a
series of macro templates which are expanded to create the wrapper
functions. Some name templates can be specified as options.  But the
overall structure of the generated code is defined by the classes and
functions in the YAML file as well as the requirements of C++ and
Fortran syntax.

Format strings contain “replacement fields” surrounded by curly braces
``{}``. Anything that is not contained in braces is considered literal
text, which is copied unchanged to the output. If you need to include
a brace character in the literal text, it can be escaped by doubling:
``{{`` and ``}}``. [Python_Format]_

Due to a feature of YAML, if a string starts with a curly brace YAML
will interpret it as a dictionary instead of as part of the
string. To avoid this behavior, strings which start with a curly brace
should be quoted::

    name : "{fmt}"

Some macros consist of blocks of code.  YAML provides a syntax for 
add multiple lines while preserving newlines::

    C_invalid_name: |
        if (! isNameValid({cpp_var})) {{
            return NULL;
        }}

Declarations may be split across several lines::

    - decl: void Sum(int len, int *values+dimension+intent(in),
                     int *result+intent(out))



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

What files are created
----------------------

Shroud will create multiple output file which must be compiled with
C++ or Fortran compilers.

One C++ file will be created for the library and one file for each C++ class.

By default, Fortran will create one file per class similar to the way
C is handled.

If one class makes use of another class in a library,
it is necessary to put all of the class
into a single file using the *F_module_per_class* option.
Each Fortran file will only contain one module to make it easier to
create makefile dependencies using pattern rules::

    %.o %.mod : %.f


How Names are created
---------------------

Shroud attempts to provide user control of names while providing
reasonable defaults.
Each name is based on the library, class, method or variable name
in the current scope.  Most names have a template which may be used
to control how the names are generated on a global scale.  Many names
may also be explicitly specified by a field.

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


How Functions are Generated
---------------------------

This section show the format templates which are used to create code.
The names in curly parens are from the format dictionary.

The C wrapper code::

    struct s_{C_type_name};
    typedef struct s_{C_type_name} {C_type_name};

C implementation::

    {C_return_type} {C_name}({C_prototype})
    {
        // c_to_cpp for the class type
        {C_const}{cpp_class} *{C_this}obj = {c_to_cpp};

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


C++ Code
^^^^^^^^

The C wrapper uses a pointer to an opaque type *C_type_name* as the 
object instance pointer.  The C++ wrapper must first cast this into
a *cpp_class* pointer.
The class's type *c_to_cpp* field is used to cast the pointer.

Next each argument uses its type *pre_call* section to convert 
the C argument into a C++ arguments. For most types this is nothing.

In addition each argument may also have a *post_call* section.

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
       namespace: one two three
    -  Class2


Local Variable
^^^^^^^^^^^^^^

*SH_* prefix on local variables.

Results are named from *fmt.C_result* or *fmt.F_result*.

Fortran option F_result.


splicers
--------


.. [Python_Format] https://docs.python.org/2/library/string.html#format-string-syntax




