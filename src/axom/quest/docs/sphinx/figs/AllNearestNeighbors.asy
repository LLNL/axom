// To turn this Asymptote source file into an image for inclusion in
// Axom's documentation,
// 1. run Asymptote:
//    asy -f png AllNearestNeighbors.asy
// 2. Optionally, use ImageMagick to convert the white background to transparent:
//    convert AllNearestNeighbors.png -transparent white AllNearestNeighbors.png

// for the search radius
import graph;

// parameters
real picscale = 40;
real dotscale = picscale / 7.0;
bool drawGraph = false;
bool drawSearch = true;

pair[] a = new pair[6];
a[0] = (-0.5, 1.8);
a[1] = (0.9, 0.8);
a[2] = (2.2, -0.3);
a[3] = (2.9, 1.9);
a[4] = (2.6, 3.2);
a[5] = (1.2, 3.8);
a = picscale*a;

pair[] b = new pair[6];
b[0] = (3.6, 1.3);
b[1] = (4.1, 0.7);
b[2] = (5.2, 0);
b[3] = (5.6, 0.9);
b[4] = (5.3, 3.4);
b[5] = (4.4, 3.1);
b = picscale*b;

pair[] c = new pair[5];
c[0] = (4.1, 4.1);
c[1] = (4.6, 4.1);
c[2] = (5.5, 4.3);
c[3] = (4.7, 5.1);
c[4] = (4.0, 4.7);
c = picscale*c;

pair[] d = new pair[6];
d[0] = (0.8, 4.7);
d[1] = (1.1, 4.4);
d[2] = (2.6, 4.0);
d[3] = (3.5, 4.4);
d[4] = (3.8, 6.2);
d[5] = (1.4, 5.9);
d = picscale*d;

// draw graph
pen graphPen = invisible;
if (drawGraph) {
  graphPen = lightgrey;
}

int[] idx = {0, 1, 2, 3, 4, 5, 6};
for (int k : idx) {
  draw(picscale*(-0.7, k)--picscale*(6.2, k), graphPen);
}
for (int k : idx) {
  draw(picscale*(k, -0.7)--picscale*(k, 6.5), graphPen);
}

// draw bodies
path pa = a[0]--a[1]--a[2]--a[3]--a[4]--a[5]--cycle;
draw(pa);
dot(a,blue+dotscale);

path pb = b[0]--b[1]--b[2]--b[3]--b[4]--b[5]--cycle;
draw(pb);
dot(b,red+dotscale);

path pc = c[0]--c[1]--c[2]--c[3]--c[4]--cycle;
draw(pc);
dot(c,orange+dotscale);

path pd = d[0]--d[1]--d[2]--d[3]--d[4]--d[5]--cycle;
draw(pd);
dot(d,green+dotscale);

// draw search radius
if (drawSearch) {
  draw (Circle(a[4], picscale), dashed);
}
