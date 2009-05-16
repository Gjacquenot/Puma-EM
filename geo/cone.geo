lc = 0.2;
lx = 1.0;
lz = 1.0;
Point(1) = {0,0,0,lc};
Point(2) = {-lx/2,0,0,lc};
Point(3) = {0,0,lz,lc};

Line(1) = {2,3};
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{1};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{2};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{5};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{8};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{11};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{14};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{17};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
  Line{20};
}
Line Loop(26) = {9,12,15,18,21,24,3,6};
Plane Surface(27) = {26};
