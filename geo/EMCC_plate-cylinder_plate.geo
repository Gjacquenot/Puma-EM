lc = 0.01;
lx = 0.14141153679245283;
offset = 2.5 * lx;

Point(1) = {-offset,0,0,lc};
Point(2) = {-offset,-lx,0,lc};
Point(3) = {-offset - lx,0,0,lc};
Point(4) = {-offset,lx,0,lc};
Point(5) = {0,-lx,0,lc};
Point(6) = {0,lx,0,lc};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Line(3) = {4,6};
Line(4) = {6,5};
Line(5) = {5,2};
Line Loop(6) = {1,2,3,4,5};
Plane Surface(7) = {6};
