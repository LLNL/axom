Options
=======




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









Functions
=========


decl

   Function declaration.
   Parsed to extract function name, type and arguments descriptions.



Annotations
===========

constructor

   Mark function as a constructor

pure

   Sets the Fortran PURE attribute.

dimension

   Sets the Fortran DIMENSION attribute.
   Pointer argument should be passed through since it is an
   array.  *value* must be *False*


value

   If true, pass-by-value; else, pass-by-reference.

intent

   IN, OUT, INOUT
   If the argument is ``const``, the default is IN.

ptr

   Argument is a pointer

reference

   Argument is a reference
