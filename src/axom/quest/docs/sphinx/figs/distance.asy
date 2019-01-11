// To turn this Asymptote source file into an image for inclusion in
// Axom's documentation,
// 1. run Asymptote:
//    asy -f png distance.asy
// 2. Optionally, use ImageMagick to convert the white background to transparent:
//    convert distance.png -transparent white distance.png

// preamble
import graph;
settings.render = 6;
size(6cm, 0);

// axes
xaxis("$x$", xmin=0, xmax=6, ticks=LeftTicks(N=6,n=2));
yaxis("$y$", ymin=0, ymax=3.8, ticks=LeftTicks);

// polyhedron
pair a = (1, 2);
pair b = (4.1, 3.2);
pair c = (5.5, 2.0);
pair d = (5.7, 0.3);
pair e = (2.8, -0.8);

path poly = a -- b -- c -- d -- e -- cycle;
draw(poly);

pair pext = (0.8, 3.2);
pair ab = b - a;
pair apext = pext - a;
pair pextp = a + dot(apext, unit(ab)) * unit(ab);

pair pint = (4.1, 1.5);
pair bc = c - b;
pair bpint = pint - b;
pair pintp = b + dot(bpint, unit(bc)) * unit(bc);

dot(pext, blue);
dot(pextp, blue);
label("ext", pext, align=N);
draw(pext--pextp, dotted, L=Label(format(length(pext - pextp)), LeftSide));

dot(pint, blue);
dot(pintp, blue);
label("int", pint, align=W);
draw(pint--pintp, dotted, L=Label(format(-1 * length(pint - pintp)), LeftSide));

