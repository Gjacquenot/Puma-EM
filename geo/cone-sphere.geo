lc = 0.149896229;
lx = 2.38567257962;
R = lx;
Point(1) = {0,0,0,lc};
Point(2) = {-R,0,0,lc};
Point(3) = {0,R,0,lc};
Point(4) = {lx,0,0,lc};
Circle(1) = {3,1,2};
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{1};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{2};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{5};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{8};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{11};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{14};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{17};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{20};
}
Line(26) = {4,3};
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{26};
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
