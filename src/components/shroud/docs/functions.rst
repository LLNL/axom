Functions
=========

What files are created
----------------------

Shroud will create multiple output file which must be compiled with
C++ or Fortran compilers.

One C++ file will be created for the library and one file for each C++ class.

By default, Fortran will create one file per class similar to the way
C is handled.
If one class makes use of another class in a library,
it is necessary to put all of the classes
into a single file using the *F_module_per_class* option.

Each Fortran file will only contain one module to make it easier to
create makefile dependencies using pattern rules::

    %.o %.mod : %.f


How Names are created
---------------------

Shroud attempts to provide user control of names while providing
reasonable defaults.
Each name is based on the library, class, method or argument name
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

    library: library
    namespace: library

    options:
      F_name_impl_template: "{library}_{underscore_name}{function_suffix}"

    function:
    -  decl: void initialize


How Functions are Wrapped
-------------------------

As each function declartion is parsed a format dictionary is created
with fields to describe the function and its arguments.
The fields are then expanded into the function wrapper.

C wrapper
^^^^^^^^^

C wrapper::

    {C_return_type} {C_name}({C_prototype})
    {
        {C_code}
    }

The wrapper is within an ``extern "C"`` block so that **C_name** will
not be mangled by the C++ compiler.

The *C_code* field has a default value of::

    {C_pre_call}
    {C_call_code}
    {C_post_call_pattern}
    {C_post_call}
    {C_return_code}

* **C_pre_call** code used to convert arguments from C to C++.
  Derived from argument types c_statements code.

* **C_call_code** code used to call the function.
  The *constructor* and *destructor* annotations will use ``new`` and ``delete``.

* **C_post_call_pattern** code from the *C_error_pattern*.
  Can be used to deal with error values.

* **C_post_call** code used with *intent(out)* arguments.

* **C_return_code** returns a value from the wrapper.


.. wrapc.py   Wrapc.write_header

C++ classes create an opaque typedef in the header file for each class::

    struct s_{C_type_name};
    typedef struct s_{C_type_name} {C_type_name};



Fortran wrapper
^^^^^^^^^^^^^^^

The template for Fortran code showing names which may 
be controlled directly by the input file::

    module {F_module_name}

      ! use_stmts
      implicit none

      type {F_derived_name}
        type(C_PTR) {F_derived_member}
      contains
        procedure :: {F_name_function} => {F_name_impl}
        generic :: {F_name_generic} => {F_name_function}, ...
      end type {F_derived_name}

      interface
        {F_C_pure_clause} {F_C_subprogram} {F_C_name}
             {F_C_result_clause} bind(C, name="{C_name}")
          ! arg_f_use
          implicit none
          ! arg_c_decl
        end {F_C_subprogram} {F_C_name}
      end interface

      interface {F_name_generic}
        module procedure {F_name_impl}
      end interface {F_name_generic}

    contains

      {F_subprogram} {F_name_impl}
        ! arg_f_use
        ! arg_f_decl
        ! pre_call
        {F_code}
        ! post_call
      end {F_subprogram} {F_name_impl}

    end module {F_module_name}


Helper functions
----------------

Functions can be created in the Fortran wrapper which have no
corresponding function in the C++ library.  This may be necessary to
add functionality which may unnecessary in C++.  For example, a
library provides a function which returns a string reference to a
name.  If only the length is desired no extra function is required in
C++ since the length is extracted used a ``std::string`` method::

    ExClass1 obj("name")
    int len = obj.getName().length();

Calling the Fortran ``getName`` wrapper will copy the string into a
Fortran array but you need the length first to make sure there is
enough room.  You can create a Fortran wrapper to get the length
without adding to the C++ library::

    classes:
      - name: ExClass1
        methods:
          - decl: int GetNameLength() const
            C_code: |
              {C_pre_call}
              return {CPP_this}->getName().length();

The generated C wrapper will use the *C_code* provided for the body::

    int AA_exclass1_get_name_length(const AA_exclass1 * self)
    {
        const ExClass1 *SH_this =
            static_cast<const ExClass1 *>(static_cast<const void *>(self));
        return SH_this->getName().length();
    }

The *C_pre_call* format string is generated by Shroud to convert the
``self`` argument into *CPP_this* and must be included in *C_code*
to get the definition.


.. Fortran shadow class


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
    } // namespace two
    } // namespace one

    class Class3 {
    };

The YAML file would look like::

    namespace: one two

    classes:
    -  Class1
       namespace: one two three
    -  Class2
    -  Class3
       namespace: -none

If a namespace starts with a ``-``, then it will be ignored.  This
allows a library to have a default namespace but have a class have no
namespace.

Local Variable
^^^^^^^^^^^^^^

*SH_* prefix on local variables.

Results are named from *fmt.C_result* or *fmt.F_result*.


Patterns
--------

To address the issue of semantic differences between Fortran and C++,
*patterns* may be used to insert additional code.  A *pattern* is a 
code template which is inserted at a specific point in the wrapper.
They are defined in the input YAML file::

   functions:
     - decl: const string& getString2+len=30()
       C_error_pattern: C_invalid_name

    patterns:
        C_invalid_name: |
            if ({cpp_var}.empty()) {{
                return NULL;
            }}

The **C_error_pattern** will insert code after the call to the C++
function in the C wrapper and before any post_call sections from the
types. The bufferified version of a function will append
``_as_buffer`` to the **C_error_pattern** value.  The *pattern* is
formated using the context of the return argument if present,
otherwise the context of the function is used.  This means that
*c_var* and *c_var_len* refer to the argument which is added to
contain the function result for the ``_as_buffer`` pattern.

The function ``getString2`` is returing a ``std::string`` referrence.
Since C and Fortran cannot deal with this directly, the empty string
is converted into a ``NULL`` pointer::
will blank fill the result::

    const char * STR_get_string2()
    {
        const std::string & SH_rv = getString2();
        // C_error_pattern
        if (SH_rv.empty()) {
            return NULL;
        }
        const char * XSH_rv = SH_rv.c_str();
        return XSH_rv;
    }



Splicers
--------

No matter how many features are added to Shroud there will always exist
cases that it does not handle.  One of the weaknesses of generated
code is that if the generated code is edited it becomes difficult to
regenerate the code and perserve the edits.  To deal with this
situation each block of generated code is surrounded by 'splicer'
comments::

    const char * STR_get_char3()
    {
    // splicer begin function.get_char3
        const char * SH_rv = getChar3();
        return SH_rv;
    // splicer end function.get_char3
    }

These comments delineate a section of code which can be replaced by
the user.  The splicer's name, ``function.get_char3`` in the example,
is used to determine where to insert the code.
In a separate file, add the begin and end splicer comments,
then add the code which should be inserted into the wrapper.  Multiple
splicer can be added to an input file.  Any text that is not within a
splicer block is ignored.  Splicers must be sorted by language.  If
the input file ends with ``.f`` or ``.f90`` it is processed as
splicers for the generated Fortran code.  Code for the C wrappers must
end with any of ``.c``, ``.h``, ``.cpp``, ``.hpp``, ``.cxx``,
``.hxx``, ``.cc``, ``.C``::

    // splicer begin function.get_char3
        const char * SH_rv = getChar3();
        SH_rv[0] = 'F';    // replace first character for Fortran
        return SH_rv + 1;
    // splicer end function.get_char3

The file with the custom splicers is added to the Shroud command line
along with the YAML file.

In addition to replacing code for a function wrapper, there are 
splicers that are generated which allow a user to insert additional
code::

    ! file_top
    module {F_module_name}
       ! module_use
       implicit none
       ! module_top

       type class1
         ! class.{cpp_class}.component_part
       contains
         ! class.{cpp_class}.generic.{F_name_generic}
         ! class.{cpp_class}.type_bound_procedure_part
       end type class1

       interface
          ! additional_interfaces
       end interface

       contains

       ! function.{F_name_function}

       ! {cpp_class}.method.{F_name_function}

       ! additional_functions

    end module {F_module_name}

.. from _create_splicer

C header::

    // class.{class_name}.CXX_declarations

    extern "C" {
    // class.{class_name}.C_declarations
    }

C implementation::

    // class.{class_name}.CXX_definitions

    extern "C" {
      // class.{class_name}.C_definitions

      // function.{underscore_name}{function_suffix}

      // class.{cpp_class}.method.{underscore_name}{function_suffix}

    }

The splicer comments can be eliminated by setting the option
**show_splicer_comments** to false. This may be useful to 
eliminate the clutter of the splicer comments.

Debugging
---------

Shroud generates a JSON file with all of the input from the YAML
and all of the format dictionaries and type maps.
This file can be useful to see which format keys are available and
how code is genenerated.

