******************************************************
1 Purpose of the guidelines and general considerations
******************************************************

These pages contain code development guidelines for the ASC Simulation 
CS Toolkit project at Lawrence Livermore National Laboratory (LLNL). 
The guidelines apply to C++ code primarily, which is the main implementation 
language for the toolkit. The majority of the guidelines herein were 
generated from the cited references, often with modifications and 
simplifications.

The guidelines emphasize code readability, correctness, portability, and 
interoperability. Agreement on coding style and following common idioms 
and patterns provides many benefits to a project with multiple developers. 
A uniform "look and feel" makes it easier to read and understand source code, 
which increases team productivity and reduces confusion and coding errors 
when developers work with code they did not write. Also, guidelines 
facilitate code reviews by enforcing consistency and focusing developers on 
common concerns. Some guidelines herein, such as naming conventions, are 
arbitrary, but all are based on practical experience and widely accepted 
sound practices. For brevity, most guidelines contain little detailed 
explanation or justification. 

The guidelines are not meant to be static, but rather to evolve with the 
needs of project developers and should be changed when needed. Changes should
be based on discussion by the development team using their collective 
professional judgment. Also, the benefits of consistency should be balanced 
with variations that allow for individual stylistic preferences and which may
be preferable in specific code circumstances. When changes are made, these 
guidelines must be updated accordingly.

Each guideline in this document is qualified with: "must", "should", or "may". 

* A "must" item is an absolute requirement. 
* A "should" item is a strong recommendation. 
* A "may" item is a potentially beneficial stylistic suggestion. 

While items having "should" and "may" qualifiers may be situation-dependent, 
they tend to enhance code readability and help reduce errors.

Two important notes:

1. The CS Toolkit contains software developed from scratch for the Toolkit 
   as well as code adopted from other sources, which was developed 
   independently of the Toolkit. These guidelines apply to software developed 
   for the Toolkit specifically. Modifying other Toolkit software to be 
   compliant with these guidelines should typically be done only if a 
   significant rewrite is undertaken for other reasons.
2. This guide is not a C++ language tutorial. Developers should be familiar 
   with the language and educate themselves on various topics as needed.
