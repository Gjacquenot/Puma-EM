lc = 0.01;
lx = 0.2;
ly = 0.1;
lz = 0.1;
Point(1) = {0,0,0,lc};
Point(2) = {0,-ly/2,0,lc};
Point(3) = {0,ly/2,0,lc};

Line(1) = {2,1};
Extrude {-lx/2,0,0} {
  Line{1};
}
Extrude {lx/2,0,0} {
  Line{1};
}
Extrude {0,0,-lz/2} {
  Line{1};
}
Extrude {0,0,lz/2} {
  Line{1};
}
Line(18) = {1,3};
Extrude {-lx/2,0,0} {
  Line{18};
}
Extrude {lx/2,0,0} {
  Line{18};
}
Extrude {0,0,-lz/2} {
  Line{18};
}
Extrude {0,0,lz/2} {
  Line{18};
}
Extrude {0,0,-lz/2} {
  Line{8};
}
Extrude {0,0,-lz/2} {
  Line{4};
}
Extrude {0,0,lz/2} {
  Line{8};
}
Extrude {0,0,lz/2} {
  Line{4};
}
