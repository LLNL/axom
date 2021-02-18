.. _inlet_verification_page_label:

############
Verification
############

Before input file data can be accessed, it must first be verified by calling the ``verify()``
method of the top-level ``Inlet`` object. This will return a ``bool`` indicating whether the
provided input conformed to the schema, and specific violations of the schema are logged via
SLIC.

This section describes the verification rules that apply to each possible element of the Inlet
hierarchy, namely, ``Table``, ``Field``, and ``Function``.

Table Verification
------------------

If a ``Table`` is marked as required (via the ``required()``) method, then the ``Table`` must be nontrivial,
meaning it either has a ``Field`` or ``Function`` that was present in the input, or contains a sub-``Table``
that does.

If a ``Table`` corresponds to an array or dictionary, the elements of the array must all be of the requested
type.  This restriction applies even if the array/dictionary was not marked as ``required``.

If a verification lambda was provided via the ``registerVerifier()`` method, this lambda must 
return ``true`` when passed the corresponding ``Table`` object.

.. note::
  Since a ``Table`` represents internal nodes in the Inlet hierarchy, its verification status is
  dependent on that of its child objects.  A ``Table`` is only considered valid if all of its child
  ``Table``, ``Field``, and ``Function`` objects are also valid.

Field Verification
------------------

If a ``Field`` is marked as required (via the ``required()``) method, then a value must be provided in the input.

Providing a value of the wrong type will result in verification failure, even if the field was not marked as ``required``.

If a range (inclusive upper/lower bounds) of valid values is specified with the ``range()`` method, both the provided value
(if applicable) and default value (if specified with ``defaultValue()``) must fall within the range.

If a set of valid values is specified with the ``validValues()`` method, both the provided value
(if applicable) and default value (if specified with ``defaultValue()``) must be included in the set.

If a verification lambda was provided via the ``registerVerifier()`` method, this lambda must 
return ``true`` when passed the corresponding ``Field`` object.

Function Verification
---------------------

If a ``Function`` is marked as required (via the ``required()``) method, then a function must be provided in the input.

If a verification lambda was provided via the ``registerVerifier()`` method, this lambda must 
return ``true`` when passed the corresponding ``Function`` object.
