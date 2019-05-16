// To turn this Asymptote source file into an image for inclusion in
// Axom's documentation, run Asymptote according to the instructions 
// further down in the file.

// preamble
size(16cm, 0);

pair[] pts = {
  (.13, .88),
  (.26, .87),
  (.11, .77),
  (.18, .78),
  (.13, .74),
  (.37, .75),
  (.12, .61),
  (.25, .51),
  (.11, .44),
  (.26, .40),
  (.12, .25),
  (.85, .38),
  (.94, .37),
  (.84, .26),
  (.92, .27),
  (.96, .28),
  (.84, .16),
  (.92, .16),
  (.93, .09)
};

// draw invisible dots to make sure we get the entire drawing
dot((-0.01, -0.01), invisible);
dot((1.01, 1.01), invisible);

// draw triangles
draw(pts[1] -- pts[0] -- pts[3] -- cycle, blue);
draw(pts[0] -- pts[2] -- pts[3] -- cycle, blue);
draw(pts[1] -- pts[3] -- pts[5] -- cycle, blue);
draw(pts[3] -- pts[2] -- pts[4] -- cycle, blue);
draw(pts[3] -- pts[4] -- pts[5] -- cycle, blue);
draw(pts[2] -- pts[6] -- pts[4] -- cycle, blue);
draw(pts[5] -- pts[4] -- pts[6] -- cycle, blue);
draw(pts[5] -- pts[6] -- pts[7] -- cycle, blue);
draw(pts[7] -- pts[6] -- pts[8] -- cycle, blue);
draw(pts[7] -- pts[8] -- pts[9] -- cycle, blue);
draw(pts[9] -- pts[8] -- pts[10] -- cycle, blue);
draw(pts[11] -- pts[13] -- pts[14] -- cycle, blue);
draw(pts[12] -- pts[11] -- pts[14] -- cycle, blue);
draw(pts[12] -- pts[14] -- pts[15] -- cycle, blue);
draw(pts[14] -- pts[13] -- pts[16] -- cycle, blue);
draw(pts[14] -- pts[16] -- pts[17] -- cycle, blue);
draw(pts[16] -- pts[18] -- pts[17] -- cycle, blue);

path[] boxes = unitsquare ^^ (0.5, 0) -- (0.5, 1) ^^ (0, 0.5) -- (1, 0.5);

// draw nth level of octree
// Set int n below, replace 0 with n and run
//    asy -f png -o showSpatialOctree0 SpatialOctree.asy
draw(unitsquare);

// draw first subdivision of octree
// Uncomment the following and run
//    asy -f png -o showSpatialOctree1 SpatialOctree.asy
real sc = 0.008;
transform shrink = shift(0.5, 0.5)*scale(1-sc)*shift(-0.5, -0.5);
draw(shrink*boxes, orange);

// draw second subdivision of octree
// Uncomment the following and run
//    asy -f png -o showSpatialOctree2 SpatialOctree.asy
sc = sc + 0.008;
shrink = shift(0.5, 0.5)*scale(1-sc)*shift(-0.5, -0.5);
draw(scale(0.5)*shift(0, 0)*shrink*boxes, blue);
draw(scale(0.5)*shift(1, 0)*shrink*boxes, blue);
draw(scale(0.5)*shift(0, 1)*shrink*boxes, blue);

// draw third subdivision of octree
// Uncomment the following and run
//    asy -f png -o showSpatialOctree3 SpatialOctree.asy
sc = sc + 0.008;
shrink = shift(0.5, 0.5)*scale(1-sc)*shift(-0.5, -0.5);
draw(scale(0.25)*shift(3, 0)*shrink*boxes, green);
draw(scale(0.25)*shift(0, 1)*shrink*boxes, green);
draw(scale(0.25)*shift(1, 1)*shrink*boxes, green);
draw(scale(0.25)*shift(3, 1)*shrink*boxes, green);
draw(scale(0.25)*shift(0, 2)*shrink*boxes, green);
draw(scale(0.25)*shift(1, 2)*shrink*boxes, green);
draw(scale(0.25)*shift(0, 3)*shrink*boxes, green);
draw(scale(0.25)*shift(1, 3)*shrink*boxes, green);
