lc = 0.1;
lx = 0.1;
ly = 0.1;
lz = 0.1;
Point(1) = {-lx/2.0,-ly/2.0,-lz/2.0,lc};
Point(2) = {lx/2.0,-ly/2.0,-lz/2.0,lc};
Point(3) = {lx/2.0,ly/2.0,-lz/2.0,lc};
Point(4) = {-lx/2.0,ly/2.0,-lz/2.0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Extrude Surface {6, {0.0,0.0,lz}};
Extrude {0,0,lz} {
  Line{11};
}
