Frequently asked questions
==========================

How do I control the names of generated functions?
--------------------------------------------------

Each C++ function may create several additional names.  Each name has a
default value describe by a option template which is inherited by
nested option blocks.  This value can be overridden for a single
function using an explicit value.

+------------------------+---------------------------------+------------------+
| Description            | Option                          | Override         |
+========================+=================================+==================+
| C wrapper              | *C_name_function_template*      | *C_name*         |
| implementation         | *C_name_method_template*        |                  |
+------------------------+---------------------------------+------------------+
| Fortran BIND(C)        | C_name.lower()                  | *F_C_name*       |
| interface              |                                 |                  |
+------------------------+---------------------------------+------------------+
| Fortran wrapper        | *F_name_impl_function_template* | *F_name_impl*    |
| implementation         | *F_name_impl_method_template*   |                  |
+------------------------+---------------------------------+------------------+
| Fortran type method    | *F_name_method_template*        | *F_name_method*  |
+------------------------+---------------------------------+------------------+
| Fortran generic name   | *F_name_generic_template*       | *F_name_generic* |
+------------------------+---------------------------------+------------------+

For *C_name* and *F_name_impl*, the template is dependent on if the function is
defined as a method of a class or not.

For example, a library has an ``initialize`` function which is
in a namespace.  In C++ it is call as::

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

