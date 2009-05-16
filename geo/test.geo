// parameters
lc = 0.0141411536792;
lx = 13.2;
ly = 13.2;
lz = 13.2;
w = 0.03;

// and now the geometry
Point(1) = {-w/2, -w/2, 0, lc};
Point(2) = {w/2, -w/2, 0, lc};
Point(3) = {w/2, w/2, 0, lc};
Point(4) = {-w/2, w/2, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Extrude {lx, 0, 0} {
  Line{2};
}
Extrude {0, ly, 0} {
  Line{3};
}
Extrude {0, 0, lz} {
  Line{4};
}
