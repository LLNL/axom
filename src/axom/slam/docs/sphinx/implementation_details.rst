******************************************************
Implementation details
******************************************************

Policy-based design 
-------------------

Handling the combinatorial explosion of features; avoid paying for what we don't need
* SizePolicy, StridePolicy, OffsetPolicy (compile time vs. runtime)
* IndirectionPolicy (none, C-array, std::vector, custom, e.g. mfem::Array)
* SubsettingPolicy (none, virtual parent, concrete parent)
* OwnershipPolicy (local, sidre, other repository)
 
(ongoing) how we simplify setup/usage for user
----------------------------------------------

* Chained initialization using named-parameter idiom
* [Generator classes to simplify types
