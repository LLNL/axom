// To turn this Asymptote source file into an image for inclusion in
// Axom's documentation,
// 1. run Asymptote:
//    asy -f png pic.asy
// 2. Optionally, use ImageMagick to convert the white background to transparent:
//    convert pic.png -transparent white pic.png

// preamble
import graph;
settings.render = 6;
size(6cm, 0);

// axes
// xaxis("$x$", xmin=0, xmax=6, ticks=LeftTicks(N=6,n=2));
// yaxis("$y$", ymin=0, ymax=3.8, ticks=LeftTicks);

// warpy polyhedron
pair a = (0.15, 1.3);
pair b = (1.75, 1.75);
pair c = (-0.5, 0.7);
pair d = (0.6, 0.9);
pair e = (1.8, 1.2);
pair f = (2.9, 1.25);
pair g = (-0.4, -0.5);
pair h = (0.6, -0.5);
pair i = (2.15, 0);
pair j = (2.9, 0.4);
pair k = (0.6, -0.7);
pair m = (2.2, -0.6);

pair da = (1, -0.8);
pair db = (0.14, -1);
pair dc = (1, -0.4);
pair dde = (1, 0.46);
pair def = (1, 0.03);
pair df = (1, 0.2);
pair dg = (1, 0);
pair dhi = (1, 0);
pair dij = (1, 0.5);
pair dj = (1, 0.1);
pair ddh = (0.8, -1);
pair dhk = (-0.1, -1);
pair dk = (0, -1);
pair dei = (0.18, -1);
pair dim = (0.19, -1);
pair dm = (0.11, -1);

path cf = c{dc}..d{dde}..e{def}..f{df};
path gj = g{dg}..h{dhi}..i{dij}..j{dj};
path ak = a{da}..d{ddh}..h{dhk}..k{dk};
path bm = b{db}..e{dei}..i{dim}..m{dm};

real blendx = 0.4;
real blendy = 0.3;

pair n = point(ak, 1 + blendy);
pair p = point(bm, 1 + blendy);
pair r = point(cf, 1 + blendx);
pair s = point(gj, 1 + blendx);

pair dn = (1-blendy)*dde + blendy*dhi;
pair dp = (1-blendy)*def + blendy*dij;
pair dr = (1-blendx)*ddh + blendx*dei;
pair ds = (1-blendx)*dhk + blendx*dim;

path np = n{dn}..p{dp};
path rs = r{dr}..s{ds};

pair o = intersectionpoint(np, rs);

draw(cf);
draw(gj);
draw(ak);
draw(bm);

dot(n);
dot(p);
dot(r);
dot(s);

draw(np, dotted);
draw(rs, dotted);

dot(o, red);
