Indexing space
^^^^^^^^^^^^^^

The UniformGrid and BVHTree spatial indexes use helper functions to address
points in :math:`R^n`.  The Morton index, or 
`Z-order curve <https://en.wikipedia.org/wiki/Z-order_curve>`_, maps all
of a region of interest in :math:`R^2` or :math:`R^3` to points on a 1D
curve, which is constructed with the goal that points that are close in
the higher-dimensional domain are close on the curve.  The 
