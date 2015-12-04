// parameters
lc = 0.02;
// inner radius
r_in = 0.049149;
// outer radius
r_out = 0.0508;
// length
ly = 0.9144;

// and now the geometry
Point(1) = {r_in, 0, 0, lc};
Point(2) = {0, 0, -r_in, lc};
Point(3) = {-r_in, 0, 0, lc};
Point(4) = {0, 0, r_in, lc};
Point(5) = {r_out, 0, 0, lc};
Point(6) = {0, 0, -r_out, lc};
Point(7) = {-r_out, 0, 0, lc};
Point(8) = {0, 0, r_out, lc};
Point(9) = {0, 0, 0, lc};


Circle(1) = {1, 9, 2};
Circle(2) = {2, 9, 3};
Circle(3) = {3, 9, 4};
Circle(4) = {4, 9, 1};
Circle(5) = {5, 9, 6};
Circle(6) = {6, 9, 7};
Circle(7) = {7, 9, 8};
Circle(8) = {8, 9, 5};
Extrude {0, -ly/2, 0} {
  Line{8, 7, 6, 5, 1, 2, 3, 4};
}
Extrude {0, ly/2, 0} {
  Line{4, 3, 2, 1, 5, 8, 7, 6};
}
Line(73) = {32, 24};
Line(74) = {33, 27};
Line(75) = {31, 28};
Line(76) = {29, 26};
Line(77) = {17, 14};
Line(78) = {15, 12};
Line(79) = {10, 23};
Line(80) = {13, 20};
Line Loop(81) = {61, 76, -41, -73};
Plane Surface(82) = {81};
Line Loop(83) = {57, 75, -53, -76};
Plane Surface(84) = {83};
Line Loop(85) = {69, 74, -49, -75};
Plane Surface(86) = {85};
Line Loop(87) = {65, 73, -45, -74};
Plane Surface(88) = {87};
Line Loop(89) = {9, -78, -37, -79};
Plane Surface(90) = {89};
Line Loop(91) = {21, -77, -25, 78};
Plane Surface(92) = {91};
Line Loop(93) = {17, 80, -29, 77};
Plane Surface(94) = {93};
Line Loop(95) = {33, -79, -13, 80};
Plane Surface(96) = {95};
