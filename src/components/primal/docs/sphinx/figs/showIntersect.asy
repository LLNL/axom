// preamble
settings.prc = false;
settings.render = 0;
import three;
size(6cm, 0);

// axes
draw(O -- 1.7X, arrow=Arrow3(DefaultHead2), L=Label("$x$", position=EndPoint, align=W));
draw(O -- 2.4Y, arrow=Arrow3(), L=Label("$y$", position=EndPoint));
draw(O -- 2Z, arrow=Arrow3(), L=Label("$z$", position=EndPoint));

// triangle 1
path3 tri1 = (1.2,0,0)--(0,1.8,0)--(0,0,1.4)--cycle;

// triangle 2
path3 tri2 = (0,0,0.5)--(0.8,0.1,1.2)--(0.8,1,1.2)--cycle;

// ray
path3 ray = (0.4,0.4,0)--(0.4,0.4,1.8);

// polygon of intersection between bbox and triangle
path3 pgon = (0.714286,0.6,0.1)--(0.457143,0.6,0.4)--(0.8,0.0857143,0.4)--(0.8,0.471429,0.1)--cycle;

// draw bounding box and other geometry
draw(box((0.1,-0.23,0.1), (0.8,0.5,0.4)), blue);
draw(pgon, deepblue);

draw(ray, arrow=Arrow3(DefaultHead2), red);
dot((0.4,0.4,0.1), red);
dot((0.4,0.4,0.622222), red);
dot((0.4,0.4,0.85), red);
draw(tri1);
draw(tri2, blue);
draw((0.420779,0.0525974,0.868182)--(0.298618,0.373272,0.76129), deepblue);
