ImplicitGrid
^^^^^^^^^^^^

Where the ``UniformGrid`` divides a rectangular region of interest into
bins, the ``ImplicitGrid`` divides each axis of the region of interest into
bins.  Each ``UniformGrid`` bin holds a list of indexes into an array,
indicating items that intersect that bin; each ``ImplicitGrid`` bin
holds a bitset indicating which item intersects that bin.

(insert figure)

The ``ImplicitGrid`` is designed for quick indexing and searching over
a static index space in a relatively coarse grid.  The following example
shows the use of the ``ImplicitGrid``:

(insert example)
