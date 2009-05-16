lc = 0.0190950610191;
angle1 = 22.62/180 * Pi;
Ssin = (Sin(angle1))^2;
cos = Cos(angle1);
length = 0.254;
half_length = length/2.0;
x1 = 4.0 * 0.0254;
x2 = 3.0 * 0.0254;
x3 = 2.0 * 0.0254;
x4 = 1.0 * 0.0254;
x5 = 0.0;

y1 = ( Sqrt(1 - (x1/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y2 = ( Sqrt(1 - (x2/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y3 = ( Sqrt(1 - (x3/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y4 = ( Sqrt(1 - (x4/half_length)^2 * Ssin ) - cos ) / (1 - cos) * 0.0254;
y5 = 1.0 * 0.0254;

Point(1) = {length/2.0,0.,0,lc/2};
Point(2) = {x1, y1, 0, lc};
Point(3) = {x2, y2, 0, lc};
Point(4) = {x3, y3, 0, lc};
Point(5) = {x4, y4, 0, lc};
Point(6) = {x5, y5, 0, lc};
CatmullRom(1) = {6,5,4,3,2,1};
Extrude Line {1, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {2, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {5, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {8, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {11, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {14, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {17, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {20, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};

// second part of the ogive
angle2 = 46.6/180 * Pi;
Ssin2 = (Sin(angle2))^2;
cos2 = Cos(angle2);
quarter_length = half_length/2.0;

x6 = -2.5 * 0.0254;
y6 = 0.0;
x7 = -2.0 * 0.0254;
y7 = ( Sqrt(1 - (x7/quarter_length)^2 * Ssin2 ) - cos2 ) / (1 - cos2) * 0.0254;
x8 = -1.5 * 0.0254;
y8 = ( Sqrt(1 - (x8/quarter_length)^2 * Ssin2 ) - cos2 ) / (1 - cos2) * 0.0254;
x9 = -1.0 * 0.0254;
y9 = ( Sqrt(1 - (x9/quarter_length)^2 * Ssin2 ) - cos2 ) / (1 - cos2) * 0.0254;
x10 = -0.5 * 0.0254;
y10 = ( Sqrt(1 - (x10/quarter_length)^2 * Ssin2 ) - cos2 ) / (1 - cos2) * 0.0254;

Point(70) = {x6, y6, 0, lc/2.0};
Point(80) = {x7, y7, 0, lc};
Point(90) = {x8, y8, 0, lc};
Point(100) = {x9, y9, 0, lc};
Point(110) = {x10, y10, 0, lc};
CatmullRom(26) = {70,80,90,100,110,6};
Extrude Line {26, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {27, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {30, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {33, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {36, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {39, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {42, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {45, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
Extrude Line {11, {1.0,0.0,0.0}, {0.0,0.0,0.0}, Pi/4};
