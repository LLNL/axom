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

=========================================================
1 Changing Existing Code
=========================================================

-----------------------------------
Follow existing code style
-----------------------------------

1.1 When modifying existing code, the style conventions already in
use in each file **must** be followed unless the scope of changes makes 
sense (see next item). This is not intended to
stifle personal creativity - mixing style is disruptive and 
may cause confusion for users and fellow developers.

1.2 When making stylistic changes to existing code, those changes **should** 
extend to a point where the style is consistent across a reasonable scope. 
This may mean that an entire file is changed to prevent multiple conflicting 
styles.

--------------------------------------------------------
Only change code from other sources when it makes sense
--------------------------------------------------------

1.3 The CS Toolkit may contain code pulled in from sources outside the 
Toolkit. These guidelines apply to code developed within the Toolkit 
primarily. The decision to modify externally-developed code that we pull 
in to the Toolkit will be evaluated on a case-by-case basis. Modifying such 
code to be compliant with these guidelines should typically be done only if 
a significant rewrite is undertaken for other reasons.
