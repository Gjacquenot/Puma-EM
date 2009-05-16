lc = 0.0999308193333;

H = 0.06;

Point(1) = {0,H,0,lc};
Point(2) = {0.03,0,-0.01,lc};
Point(3) = {-0.03,0,0.01,lc};
Point(4) = {-0.02,0,-0.04,lc};

Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,4};
Line(4) = {2,3};
Line(5) = {3,4};
Line(6) = {4,2};

Line Loop(1) = {2,-4,-1};
Plane Surface(1) = {1};
Line Loop(2) = {-2,3,-5};
Plane Surface(2) = {2};
Line Loop(3) = {1,-6,-3};
Plane Surface(3) = {3};
Line Loop(4) = {4,5,6};
Plane Surface(4) = {4};

Physical Surface(1)={1,2,3,4}; 
