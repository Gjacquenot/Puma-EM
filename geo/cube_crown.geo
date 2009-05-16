lc = 0.1;
lx = 0.1;
w = 0.02; // width of crown
Point(1) = {-lx/2.0,-lx/2.0,0.0,lc};
Point(2) = {lx/2.0,-lx/2.0,0.0,lc};
Point(3) = {lx/2.0,lx/2.0,0.0,lc};
Point(4) = {-lx/2.0,lx/2.0,0.0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {2,3,4,1};


Point(5) = {-(lx+w)/2.0,-(lx+w)/2.0,0.0,lc};
Point(6) = {(lx+w)/2.0,-(lx+w)/2.0,0.0,lc};
Point(7) = {(lx+w)/2.0,(lx+w)/2.0,0.0,lc};
Point(8) = {-(lx+w)/2.0,(lx+w)/2.0,0.0,lc};

Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};
Extrude {0, 0, lx/2.0} {
  Line{3, 4, 1, 2};
}
Extrude {0, 0, -lx/2.0} {
  Line{1, 2, 3, 4};
}
Line Loop(42) = {14, 18, 22, 10};
Plane Surface(43) = {42};
Line Loop(44) = {26, 30, 34, 38};
Plane Surface(45) = {44};
Line Loop(46) = {6, 7, 8, 9};
Plane Surface(47) = {46, 5};
