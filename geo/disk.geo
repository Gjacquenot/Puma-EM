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


Line(5) = {1, 4};
Line(6) = {1, 3};
Line(7) = {1, 2};
Line(8) = {1, 5};
Line Loop(9) = {5, -2, -6};
Plane Surface(10) = {9};
Line Loop(11) = {6, -1, -7};
Plane Surface(12) = {11};
Line Loop(13) = {7, -4, -8};
Plane Surface(14) = {13};
Line Loop(15) = {8, -3, -5};
Plane Surface(16) = {15};
