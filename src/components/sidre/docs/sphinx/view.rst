.. ##
.. ## Copyright (c) 2017-18, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

.. _view-label:

==========
View
==========

A Sidre view describes data and provides access to it. 

.. note:: View objects can only be created and destroyed using Group methods 
          provided for this. The View constructor and destructor are private.

Each View has a name and is owned by one Group in a Sidre Group hierarchy;
its *owning* Group. A View maintains a pointer to the Group that owns it.

.. note:: * The name (string) of a View **must be unique** within its
            owning Group.
          * A View has a unique integer identifier within its owning group, 
            which is generated when the View is created.
          * Views in a Group can be accessed by name or integer id.

A View object can describe and provide access to data associated with a 
pointer in one of four ways described below. In that case, a View data 
description includes: a data type, a length (number of elements), an offset 
and a stride (based on the pointer address and data type). 

  * A View can describe (a subset of) data owned by an existing Buffer. 
    In this case, the Buffer is manually *attached* to the View and the
    View's data description is applied to the Buffer data. Data can be 
    (re)allocated or deallocated by the View if and only if it is the only 
    View attached to the buffer. **In general, a Buffer can be attached
    to more than one View.**
  * A View description is used to allocate data for the View using semantics 
    similar to Buffer data description and allocation (see :ref:`buffer-label`).
    In this case, the View is exclusively associated with a Buffer and no 
    other View is allowed to (re)allocate or deallocate the data held by the 
    Buffer.
  * A View can **describe** data associated with a pointer to an *external* 
    data allocation. In this case, the View cannot (re)allocate or deallocate 
    the data. However, all other View operations can be applied to the data
    in essentially the same ways as the previous two cases.
  * A View can hold a pointer to an undescribed (*opaque*) data pointer. In 
    this case, the View knows nothing about the type or structure of the data; 
    it can only provide access to it. A user is entirely responsible for 
    casting the pointer to a proper type, knowing the length of the data, etc.

A View may also refer to a scalar quantity or a string.

Before we summarize the Sidre View interface, we present some View concepts
that describe various *states* a View can be in at any given time. Hopefully,
this will help provide some useful context for the method descriptions that 
follow.

.. figure:: figs/sidre-states.png

   This table provides a summary of the Sidre View state matrix. Each row 
   corresponds to *data association* and each column refers to a *data state*.
   The True/False entries in the cells refer to return values of the
   View method at the top of each column. The circumstances under which those
   values are returned are noted as well.
