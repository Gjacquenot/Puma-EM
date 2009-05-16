lc = 0.01;
lx = 0.1; // supposed to be the great axis
factor = 0.5; // the multiplication factor for the left side of the ellipsoid
ly = 0.02; // the target becomes a circular "cushion" if ly>lx
Point(1) = {0,0,0,lc};
Point(2) = {0,ly,0,lc};
Point(3) = {lx,0,0,lc};
Point(4) = {-lx * factor,0,0,lc};
Ellipse(1) = {2,1,3,3};
Ellipse(2) = {2,1,4,4};
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{1};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{3};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{6};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{9};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{12};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{15};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{18};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{21};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{2};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{27};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{30};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{33};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{36};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{39};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{42};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{45};
}
