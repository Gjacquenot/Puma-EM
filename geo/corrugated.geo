lx = 1.0; // lx = lambda
lc = lx/10;
totalHeight = 0.3 * lx;
grooveHeight = 0.25 * lx;
grooveWidth = 0.05 * lx;
grooveSpacing = 0.01 * lx;
totalLength = 17 * (grooveWidth + grooveSpacing) + grooveSpacing;
Point(1) = {-totalLength/2.0,-lx/2.0,grooveHeight-totalHeight,lc};
Point(2) = {totalLength/2.0,-lx/2.0,grooveHeight-totalHeight,lc};
Point(3) = {totalLength/2.0,lx/2.0,grooveHeight-totalHeight,lc};
Point(4) = {-totalLength/2.0,lx/2.0,grooveHeight-totalHeight,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};



// first corrugation
basePointNumber = 10;
baseLineNumber = 10;
basePosition = -totalLength/2.0;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};


Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// second corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// third corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// fourth corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// fifth corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// sixth corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// 7th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// 8th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// 9th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// 10th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// 11th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// 12th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// 13th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// 14th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// 15th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// 16th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// 17th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};

// 18th corrugation
basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber + 25;
basePosition = basePosition + grooveWidth + grooveSpacing;
Point(basePointNumber) = {basePosition,-lx/2.0,0,lc};
Point(basePointNumber+1) = {basePosition,lx/2.0,0,lc};
Point(basePointNumber+2) = {basePosition,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+3) = {basePosition,lx/2.0,grooveHeight,lc};
Point(basePointNumber+4) = {basePosition + grooveSpacing,-lx/2.0,grooveHeight,lc};
Point(basePointNumber+5) = {basePosition + grooveSpacing,lx/2.0,grooveHeight,lc};
Point(basePointNumber+6) = {basePosition + grooveSpacing,-lx/2.0,0,lc};
Point(basePointNumber+7) = {basePosition + grooveSpacing,lx/2.0,0,lc};
Line(baseLineNumber) = {basePointNumber,basePointNumber+1};
Line(baseLineNumber+1) = {basePointNumber+2,basePointNumber+3};
Line(baseLineNumber+2) = {basePointNumber+4,basePointNumber+5};
Line(baseLineNumber+3) = {basePointNumber+6,basePointNumber+7};
Line(baseLineNumber+4) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+5) = {basePointNumber+1,basePointNumber+3};
Line(baseLineNumber+6) = {basePointNumber+2,basePointNumber+4};
Line(baseLineNumber+7) = {basePointNumber+3,basePointNumber+5};
Line(baseLineNumber+8) = {basePointNumber+4,basePointNumber+6};
Line(baseLineNumber+9) = {basePointNumber+5,basePointNumber+7};
Line(baseLineNumber+10) = {basePointNumber,basePointNumber+6};
Line(baseLineNumber+11) = {basePointNumber+1,basePointNumber+7};

Line Loop(baseLineNumber+12) = {baseLineNumber+4,baseLineNumber+1,-(baseLineNumber+5),-baseLineNumber};
Plane Surface(baseLineNumber+13) = {baseLineNumber+12};
Line Loop(baseLineNumber+14) = {baseLineNumber+4,baseLineNumber+6,baseLineNumber+8,-(baseLineNumber+10)};
Plane Surface(baseLineNumber+15) = {baseLineNumber+14};
Line Loop(baseLineNumber+16) = {baseLineNumber+1,baseLineNumber+7,-(baseLineNumber+2),-(baseLineNumber+6)};
Plane Surface(baseLineNumber+17) = {baseLineNumber+16};
Line Loop(baseLineNumber+18) = {baseLineNumber+5,baseLineNumber+7,baseLineNumber+9,-(baseLineNumber+11)};
Plane Surface(baseLineNumber+19) = {baseLineNumber+18};
Line Loop(baseLineNumber+20) = {baseLineNumber+2,baseLineNumber+9,-(baseLineNumber+3),-(baseLineNumber+8)};
Plane Surface(baseLineNumber+21) = {baseLineNumber+20};


// now the surfaces between the dents
basePointNumber = 16;
baseLineNumber = 457;
baseLineNumber2 = 13;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

basePointNumber = basePointNumber + 8;
baseLineNumber = baseLineNumber+4;
baseLineNumber2 = baseLineNumber2+25;
Line(baseLineNumber) = {basePointNumber,basePointNumber+2};
Line(baseLineNumber+1) = {basePointNumber+1,basePointNumber+3};
Line Loop(baseLineNumber+2) = {baseLineNumber2,baseLineNumber+1,-(baseLineNumber2+22),-baseLineNumber};
Plane Surface(baseLineNumber+3) = {baseLineNumber+2};

// last elements...
Line(525) = {1,10};
Line(526) = {2,152};
Line(527) = {3,153};
Line(528) = {4,11};
Line Loop(529) = {4,525,10,-528};
Plane Surface(530) = {529};
Line Loop(531) = {1,526,-445,-521,-420,-517,-395,-513,-370,-509,-345,-505,-320,-501,-295,-497,-270,-493,-245,-489,-220,-485,-195,-481,-170,-477,-145,-473,-120,-469,-95,-465,-70,-461,-45,-457,-20,-525};
Plane Surface(532) = {531};
Line Loop(533) = {2,527,-438,-526};
Plane Surface(534) = {533};
Line Loop(535) = {3,528,21,458,46,462,71,466,96,470,121,474,146,478,171,482,196,486,221,490,246,494,271,498,296,502,321,506,346,510,371,514,396,518,421,522,446,-527};
Plane Surface(536) = {535};
