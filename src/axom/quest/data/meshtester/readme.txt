
This directory contains data files for verifying 
axom::quest::findTriMeshIntersections() and related functions.  Each test is
specified by a file with the extension ".test".  
These text files have six lines:

1. The test name, a distinctive string.
2. The STL mesh file, relative to this directory.  Suggested locations are in
   this directory or its parent.
3. An integer corresponding to the underlying values of quest::WatertightStatus
   0: the mesh is watertight (after welding vertices)
   1: the mesh is not watertight (has at least one boundary edge)
   2: the mesh has some other error, possibly an edge incident in more than
      two faces
4. An integer corresponding to the genus of the mesh after vertices have been welded
5. A list of pairs of integers, specifying the indices of intersecting
   triangles.  Each index should be the 0-based index of the triangle in the
   STL file.
6. A list of integers, specifying the STL file indices of degenerate triangles.

The integers on lines 5 and 6 should be separated by whitespace.
