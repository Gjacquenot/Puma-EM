lc = 0.0164721130769;
lx = 1.0;
lz = 1.0;

Angle = Pi/2;

Point(1) = {0,0,-lz/2,lc};
Point(2) = {Cos(Angle/2) * lx,-Sin(Angle/2) * lx,-lz/2,lc};
Point(3) = {Cos(Angle/2) * lx,Sin(Angle/2) * lx,-lz/2,lc};
Point(4) = {0,0,lz/2,lc};
Point(5) = {Cos(Angle/2) * lx,-Sin(Angle/2) * lx,lz/2,lc};
Point(6) = {Cos(Angle/2) * lx,Sin(Angle/2) * lx,lz/2,lc};


Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {4,5};
Line(4) = {4,6};
Line(5) = {6,3};
Line(6) = {4,1};
Line(7) = {5,2};
Line Loop(8) = {3,7,-1,-6};
Plane Surface(9) = {8};
Line Loop(10) = {6,2,-5,-4};
Plane Surface(11) = {10};
