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
it is necessary to put all of the class
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
with fields to describe the function.  The fields are then expanded
into the function wrapper.

C wrapper::

    {C_return_type} {C_name}({C_prototype})
    {
        {C_code}
    }

The wrapper is within an ``extern "C"`` block so that **C_name** will
not be mangled by the C++ compiler.

**C_code** can be specified in the YAML file create a function for Fortran which
has no corresponding function in the C++ library::

    - decl: int GetNameLength() const
      doxygen:
        brief: helper function for Fortran to get length of name.
      C_code: |
        {C_pre_call}
        return {CPP_this}->getName().length();




        {rv_decl} = {CPP_this_call}{method_name}{CPP_template}({C_call_list});
        // pre_call
        // post-call
        // return_line



The template for Fortran code showing names which may 
be controlled directly by the input file::

    module {F_module_name}

      type {F_derived_name}
        type(C_PTR) {F_derived_member}
      contains
        procedure :: {F_name_function} => {F_name_impl}
        generic :: {F_name_generic} => {F_name_function}, ...
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


C++ Classes
^^^^^^^^^^^

The C wrapper uses a pointer to an opaque type *C_type_name* as the 
object instance pointer.


.. wrapc.py   Wrapc.write_header

The C wrapper header file::

    struct s_{C_type_name};
    typedef struct s_{C_type_name} {C_type_name};



The C++ wrapper must first cast this into
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

Fortran option F_result.

Patterns
--------

blah blah blah


Splicers
--------

No matter how many features are added to Shroud there will always exist
cases that it does not handle.  Once of the weaknesses of generated
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
the user.  In a separate file, add the begin and end splicer comments,
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


In addition to replacing code for a function wrapper, there are many
empty splicers that are generated which allow a user to insert additional
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

      // {cpp_class}.method.{underscore_name}{function_suffix}

    }

The splicer comments can be eliminated by setting the option
**show_splicer_comments** to false.
