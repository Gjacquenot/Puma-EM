// parameters
lc = 0.02398339664;
lx = 0.151;
lz = 0.051;

// and now the geometry
Point(1) = {0, 0, 0, lc};
Point(2) = {-lx, 0, 0, lc};
Point(3) = {0, -lx, 0, lc};
Point(4) = {lx, 0, 0, lc};
Point(5) = {0, lx, 0, lc};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Extrude {0,0,lz/2.0} {
  Line{2,1,4,3};
}
Extrude {0,0,-lz/2.0} {
  Line{1,4,3,2};
}
Line Loop(37) = {9,5,17,13};
Plane Surface(38) = {37};
Line Loop(39) = {21,33,29,25};
Plane Surface(40) = {39};
