.. ##
.. ## Copyright (c) 2016, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## All rights reserved.
.. ##
.. ## This file cannot be distributed without permission and
.. ## further review from Lawrence Livermore National Laboratory.
.. ##


=====================================
7 Code Formatting
=====================================

--------------------------------------------------------------------
7.1 Conditional statements and loops
--------------------------------------------------------------------

7.1.1 Curly braces **should** be used in all conditionals, loops, etc. 
even when the content inside the braces is a "one-liner". 

       This helps prevent coding errors and misinterpretation of intent. 
       For example, this::

          if (done) { ... }

       is preferable to this::

          if (done) ...

7.1.2 One-liners **may** be used for "if" conditionals with 
"else/else if"  clauses when the resulting code is clear. 

       For example, either of the following styles **may** be used::

          if (done) {
             id = 3;
          } else {
             id = 0;
          }

       or::

          if (done) { id = 3; } else { id = 0; }

7.1.3 Complex "if/else if" conditionals with many "else if" clauses 
**should** be avoided.

      Such statements can always be refactored using local boolean variables 
      or "switch" statements. Doing so often makes code easier to read and 
      understand and may improve performance.

7.1.4 An explicit test for zero/nonzero **must** be used in a conditional 
unless the tested quantity is a boolean or pointer type. 

      For example, a conditional based on an integer value should use::

         if (num_lines != 0) {

      not::

         if (num_lines) {


--------------------------------------------------------------------
7.2 White Space and Code Alignment
--------------------------------------------------------------------

Most conventions for indentation, spacing and code alignment 
preferred by the team are enforced by using the `uncrustify` tool. 
It can be run from the top-level CS Toolkit directory...

.. note :: Show how to run uncrustify on the code and where the format
           options are defined.

Not all preferred formatting conventions are supported by uncrustify.
The following guidelines provide additional recommendations to make
code easier to read and understand.


Use White Space To Make Code Easier to Read
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.1 Blank lines and indentation **should** be used throughout code to 
enhance readability. 

      Examples of helpful white space include:

         * Between operands and operators in arithmetic expressions.
         * After reserved words, such as "while", "for", "if", "switch", etc. 
           and before the parenthesis or curly brace that follows.
         * After commas separating arguments in functions.
         * After semi-colons in for-loop expressions.
         * Before and after curly braces in almost all cases.

7.2.2 White space **must not** appear between a function name and the opening 
parenthesis to the argument list. In particular, if a function call is broken 
across source lines, the break **must not** come between the function name and 
the opening parenthesis.

7.2.3 Tabs **must not** be used for indentation since this can be problematic 
for developers with different text editor settings.


Align Code Vertically to Show Scope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.4 When function arguments (in either a declaration or implementation)
appear on multiple lines, the arguments **should** be aligned vertically 
for readability.

7.2.5 All statements within a function body **should** be indented within the surrounding curly braces.

7.2.6 All source lines in the same scope **should** be aligned vertically.
Continuation of previous lines **may** be indented if it make the code easier
to read.


Break Lines Where It Makes Sense
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.7 When a line is broken at a comma or semi-colon, it **must** be broken 
after the comma or semi-colon, not before. 

7.2.8 When a source line is broken at an arithmetic operator 
(i.e., , -, etc.), it **should** be broken after the operator, not before. 


Use Parentheses For Clarity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.9 Parentheses **should** be used in non-trivial mathematical and logical 
expressions to clearly indicate structure and intended order of operations. 
Do not assume everyone who reads the code knows all the rules for operator 
precedence.
