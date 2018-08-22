// To turn this Asymptote source file into an image for inclusion in
// Axom's documentation, run Asymptote:
//    asy -f png showUniformGrid.asy

// preamble
size(10cm, 0);

// points
pair[] p = new pair[8];
p[0] = (0.3,  0.93);
p[1] = (0.1,  0.85);
p[2] = (0.3,  0.78);
p[3] = (0.18, 0.36);
p[4] = (0.8,  0.58);
p[5] = (0.6,  0.5);
p[6] = (0.55, 0.42);
p[7] = (0.61, 0.1);

// triangles
path[] t = new path[6];
t[0] = p[0]--p[1]--p[2]--cycle;
t[1] = p[2]--p[1]--p[3]--cycle;
t[2] = p[2]--p[3]--p[6]--cycle;
t[3] = p[6]--p[3]--p[7]--cycle;
t[4] = p[4]--p[2]--p[6]--cycle;
t[5] = p[4]--p[5]--p[7]--cycle;

// invisible dots at the corners
dot((-.005, -.005), invisible);
dot((1.005, 1.005), invisible);

// uniform grid
real scfactor = 0.333;
path unitsq = (0, 0)--(1, 0)--(1, 1)--(0, 1)--cycle;
transform onethird = scale(scfactor);
for (int s = 0; s < 3; ++s) {
  for (int t = 0; t < 3; ++t) {
    transform sh = shift(scfactor * s, scfactor * t);
    draw(sh * onethird * unitsq, grey);
  }
}

// draw triangles
draw(t[0], blue);
draw(t[1], blue);
draw(t[2]);
draw(t[3]);
draw(t[4]);
draw(t[5], orange);
