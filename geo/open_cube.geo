lc = 0.0155397293179;
lx = 0.1;
Point(1) = {-lx/2.0,-lx/2.0,-lx/2.0,lc};
Point(2) = {lx/2.0,-lx/2.0,-lx/2.0,lc};
Point(3) = {lx/2.0,lx/2.0,-lx/2.0,lc};
Point(4) = {-lx/2.0,lx/2.0,-lx/2.0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Extrude {0,0,lx} {
  Line{3,4,1,2};
}
Line Loop(21) = {9,13,17,5};
Plane Surface(22) = {21};
