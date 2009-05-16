lc = 0.01;
lx = 0.1; // supposed to be the great axis
factor = 0.5; // the multiplication factor for the -x side of the elliptic cylinder
ly = 0.03;
lz = 0.2;
Point(1) = {0,0,0,lc};
Point(2) = {0,ly,0,lc};
Point(3) = {lx,0,0,lc};
Point(4) = {-lx * factor,0,0,lc};
Point(5) = {0,-ly,0,lc};
Ellipse(1) = {2,1,3,3};
Ellipse(2) = {5,1,3,3};
Ellipse(3) = {5,1,4,4};
Ellipse(4) = {2,1,4,4};
Line Loop(5) = {1,-2,3,-4};
Plane Surface(6) = {5};
Extrude {0,0,lz/2.0} {
  Surface{6};
}
Extrude {0,0,-lz/2.0} {
  Surface{6};
}
