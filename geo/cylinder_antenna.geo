// parameters
lc = 0.05;
lx = 0.051;
lz = 0.3;
frill_width = 0.02;

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
  Line{1,2,3,4};
}
Extrude {0,0,-lz/2.0} {
  Line{1,2,3,4};
}
Line Loop(37) = {5,9,13,17};
Plane Surface(38) = {37};
Line Loop(39) = {25,29,33,21};
Plane Surface(40) = {39};


Point(newp) = {-(lx+frill_width), 0, 0, lc};
Point(newp) = {0, -(lx+frill_width), 0, lc};
Point(newp) = {lx+frill_width, 0, 0, lc};
Point(newp) = {0, lx+frill_width, 0, lc};

Circle(41) = {24,1,25};
Circle(42) = {25,1,26};
Circle(43) = {26,1,27};
Circle(44) = {27,1,24};
Line(45) = {3,25};
Line(46) = {4,26};
Line(47) = {5,27};
Line(48) = {2,24};
Line Loop(49) = {46,43,-47,-3};
Plane Surface(50) = {49};
Line Loop(51) = {47,44,-48,-4};
Plane Surface(52) = {51};
Line Loop(53) = {48,41,-45,-1};
Plane Surface(54) = {53};
Line Loop(55) = {45,42,-46,-2};
Plane Surface(56) = {55};
Physical Surface("Radiating structure") = {38,20,8,12,16,28,24,36,32,40};
Physical Surface("Magnetic frill") = {56,50,52,54};
