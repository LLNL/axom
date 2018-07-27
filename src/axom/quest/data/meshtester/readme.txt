
This directory contains data files for verifying 
axom::quest::findTriMeshIntersections() and related functions.  Each test is
specified by a file with the extension ".test".  These text files have
four lines:

1. The test name, a distinctive string.
2. The STL mesh file, relative to this directory.  Suggested locations are in
   this directory or its parent.
3. A list of pairs of integers, specifying the indices of intersecting
   triangles.  Each index should be the 0-based index of the triangle in the
   STL file.
4. A list of integers, specifying the STL file indices of degenerate triangles.

The integers on lines 3 and 4 should be separated by whitespace.
