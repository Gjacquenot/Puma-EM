lc = 0.0254061405085;
angle = 22.62/180 * Pi;
Ssin = (Sin(angle))^2;
cos = Cos(angle);
length = 0.254;
half_length = length/2.0;
x1 = -4.0 * 0.0254;
x2 = -3.0 * 0.0254;
x3 = -2.0 * 0.0254;
x4 = -1.0 * 0.0254;
x5 = 0.0;

y1 = ( Sqrt(1 - (x1/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y2 = ( Sqrt(1 - (x2/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y3 = ( Sqrt(1 - (x3/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y4 = ( Sqrt(1 - (x4/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y5 = 1.0 * 0.0254;

Point(1) = {-length/2.0,0.,0,lc/4};
Point(2) = {x1, y1, 0, lc};
Point(3) = {x2, y2, 0, lc};
Point(4) = {x3, y3, 0, lc};
Point(5) = {x4, y4, 0, lc};
Point(6) = {x5, y5, 0, lc};
Point(7) = {length/2.0,0.,0,lc/4};
Point(8) = {-x1, y1, 0, lc};
Point(9) = {-x2, y2, 0, lc};
Point(10) = {-x3, y3, 0, lc};
Point(11) = {-x4, y4, 0, lc};



CatmullRom(1) = {1,2,3,4,5,6};
CatmullRom(2) = {7,8,9,10,11,6};
Extrude Line {1, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {3, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {6, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {9, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {12, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {15, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {18, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {21, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {2, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {27, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {30, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {33, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {36, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {39, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {42, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {45, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};

