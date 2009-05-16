lc = 0.0141411536792;
lx = 0.530293262972;
ly = 0.353528841981;
lz = 0.353528841981;
Point(1) = {0,-ly/2.0,-lz/2.0,lc};
Point(2) = {0,ly/2.0,-lz/2.0,lc};
Point(3) = {0,ly/2.0,lz/2.0,lc};
Point(4) = {0,-ly/2.0,lz/2.0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Extrude {-lx,0,0} {
  Line{2,3,4,1};
}
Line Loop(21) = {13,17,5,9};
Plane Surface(22) = {21};
