lc = 0.00942743578616;
lx = 0.07;
x_center = 0;
y_center = 0;
z_center = 0;

Point(1) = {x_center, y_center, z_center, lc};
Point(2) = {x_center - lx, y_center, z_center, lc};

Point(3) = {x_center,y_center-lx,z_center,lc};
Point(4) = {x_center-lx*Sqrt(2)/2,y_center-lx*Sqrt(2)/2,z_center,lc};
Point(5) = {x_center+lx*Sqrt(2)/2,y_center-lx * Sqrt(2)/2,z_center,lc};
Point(6) = {x_center+lx,y_center,z_center,lc};

// the basis circles
Circle(1) = {2,1,4};
Circle(2) = {4,1,3};


Circle(3) = {3,1,5};
Circle(4) = {5,1,6};
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{1};
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
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{23};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{2};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{29};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{33};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{37};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{41};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{45};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{49};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{53};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{3};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{61};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{65};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{69};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{73};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{77};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{81};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{85};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{4};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{93};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{96};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{99};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{102};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{105};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{108};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{111};
}
