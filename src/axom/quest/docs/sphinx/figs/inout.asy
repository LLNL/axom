// To turn this Asymptote source file into an image for inclusion in
// Axom's documentation,
// 1. run Asymptote:
//    asy -f png inout.asy
// 2. Optionally, use ImageMagick to convert the white background to transparent:
//    convert inout.png -transparent white inout.png

// preamble
settings.render = 6;
import three;
size(6cm, 0);

// axes
draw(O -- 2.1X, arrow=Arrow3(DefaultHead2), L=Label("$x$", position=EndPoint, align=N));
draw(O -- 1.7Y, arrow=Arrow3(), L=Label("$y$", position=EndPoint));
draw(O -- 1.35Z, arrow=Arrow3(), L=Label("$z$", position=EndPoint, align=W));

// tetrahedron for simplicity
triple a = (1.9, 0.2, 0);
triple b = (1.8, 1.5, 0);
triple c = (0.8, 1.1, 0);
triple d = (1.1, 0.6, 1.2);

draw(a--b--c--a--d--b);
draw(d--c);

// probe points
triple in = (1.6, 1, 0.3);
triple on = 0.65*d + 0.35*b;
triple out = (1.05, 1.45, 0.2);

dot(in, blue);
label("in", in, align=W);
draw (in--(in.x, in.y, 0), dotted);
dot(on, blue);
label("on", on, align=W);
dot(out, blue);
label("out", out, align=E);
draw (out--(out.x, out.y, 0), dotted);

