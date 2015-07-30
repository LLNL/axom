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

  Function with return a ``char *`` will instead by converted to subroutine
  which require an additional argument for the result.

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

   Fortran PURE

