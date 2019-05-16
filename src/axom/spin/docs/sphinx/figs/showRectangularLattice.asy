// To turn this Asymptote source file into an image for inclusion in
// Axom's documentation, run Asymptote:
//    asy -f png showRectangularLattice.asy

// preamble
size(10cm, 0);

real xmin = -0.8;
real xmax = 3.4;
real ymin = -1.4;
real ymax = 1.8;
real ticsize = 0.1;

path vertbar = (0, ymin) -- (0, ymax);
path horzbar = (xmin, 0) -- (xmax, 0);

// axis
draw(horzbar);
draw(vertbar);

// tic marks
path xtic = (0, 0) -- (0, ticsize);
path ytic = (0, 0) -- (ticsize, 0); 

draw(shift(1, 0)*xtic);
draw(shift(2, 0)*xtic);
draw(shift(0, -1)*ytic);
draw(shift(0, 1)*ytic);
label("0", (0, 0), SW);
label("1", (1, 0), S);
label("2", (2, 0), S);
label("-1", (0, -1), W);
label("1", (0, 1), W);

// RectangularLattice
pair origin = (-0.6, -0.2);
pair spacing = (1.2, 0.8);

// Shade the bin
path thebin = shift(origin)*scale(spacing.x, spacing.y)*shift(2, 1)*unitsquare;
fill(thebin, lightgrey);

// Draw the lines
for (int s = 0; s < 4; ++s) {
  draw(shift(0, origin.y + spacing.y*(s-1))*horzbar, orange);
  draw(shift(origin.x + spacing.x*s, 0)*vertbar, orange);
}

// Circle the origin
draw(circle(origin, 0.11), orange);

// Query points
pair q1 = (2.0, 1.2);
dot(q1);
label("$A$", q1, align=NE);

pair q2 = (2.3, 0.8);
dot(q2);
label("$B$", q2, align=E);

