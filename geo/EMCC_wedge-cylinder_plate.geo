lc = 0.0155397293179;
lx = 0.141411536792;
offset = 2.0 * lx * Cos(30.0/180.0 * Pi);

Point(1) = {-offset,0,0,lc};
Point(2) = {-offset,-lx,0,lc};
Point(3) = {-offset-lx,0,0,lc};
Point(4) = {-offset,lx,0,lc};
Point(5) = {0,0,0,lc};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Line(3) = {5,4};
Line(4) = {5,2};
Line Loop(5) = {3,-2,-1,-4};
Plane Surface(6) = {5};
