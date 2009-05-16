// parameters
lc = 0.0251926435294;
d = 9.936 * 0.0254;

/************
 * Functions
 ************/
Function QuarterEllipse1
  x = d*t * t_fact;
  y = 4.83345 * d * ( Sqrt(1 - (t*t_fact/2.08335)^2) - 0.96);
  z = 1.61115 * d * ( Sqrt(1 - (t*t_fact/2.08335)^2) - 0.96);
  Point(point_number) = {x, 0, 0, lc*lc_fact}; // ellipse center
  CenterNumber = point_number;
  Psi_array[] = {0.0, 15.0, 45.0, 90.0};
  For i In {0:3}
    point_number = newp;
    thePointNumber[i] = point_number;
    psi = Psi_array[i]/180.0*Pi ;
    Point(point_number) = {x, y*Cos(psi), z*Sin(psi), lc*lc_fact}; point_number = newp;
  EndFor
  For i In {0:2}
    Ellipse(newreg) = {thePointNumber[i],CenterNumber,thePointNumber[3],thePointNumber[i+1]};
  EndFor
Return

Function QuarterEllipse2
  x = d*t * t_fact;
  y = yfact * d * ( Sqrt(1 - (t*t_fact*3.0/1.25)^2) );
  z = zfact * d * ( Sqrt(1 - (t*t_fact*3.0/1.25)^2) );
  Point(point_number) = {x, 0, 0, lc*lc_fact}; // ellipse center
  CenterNumber = point_number;
  Psi_array[] = {0.0, 15.0, 45.0, 90.0};
  For i In {0:3}
    point_number = newp;
    thePointNumber[i] = point_number;
    psi = Psi_array[i]/180.0*Pi ;
    Point(point_number) = {x, y*Cos(psi), z*Sin(psi), lc*lc_fact}; point_number = newp;
  EndFor
  For i In {0:2}
    Ellipse(newreg) = {thePointNumber[i],CenterNumber,thePointNumber[3],thePointNumber[i+1]};
  EndFor
Return

/************
 * Geometry
 ************/
// first part of the almond
t = 1.75/3.0;

// tip of the almond
Point(1) = {d*t, 0.0, 0.0, lc/2.0};
point_number = 2;

// 0.975
lc_fact = 0.5;
t_fact = 0.975;
Call QuarterEllipse1 ;

// 0.95
lc_fact = 0.75;
t_fact = 0.95;
Call QuarterEllipse1 ;

lc_fact = 1.0;
t_fact = 0.9;
For j In {1:10}
  Call QuarterEllipse1 ;
  t_fact -= 0.1;
EndFor


// second part of the almond

t = -1.25/3.0;
yfact = 0.58/3.0;
zfact = 0.58/9.0;

lc_fact = 1.0;
t_fact = 0.1;
For j In {1:9}
  Call QuarterEllipse2 ;
  t_fact += 0.1;
EndFor


// 0.95
lc_fact = 1.0;
t_fact = 0.95;
Call QuarterEllipse2 ;

// 0.99
lc_fact = .75;
t_fact = 0.99;
Call QuarterEllipse2 ;

// tip of the almond
Point(point_number) = {d*t, 0.0, 0.0, lc/2.0};

// curves

CatmullRom(70) = {1,6,11};
CatmullRom(71) = {1,5,10};
CatmullRom(72) = {1,3,8};
CatmullRom(73) = {11,16,21};
CatmullRom(74) = {10,15,20};
CatmullRom(75) = {9,14,19};
CatmullRom(76) = {8,13,18};
CatmullRom(77) = {21,26,31,36};
CatmullRom(78) = {20,25,30,35};
CatmullRom(79) = {19,24,29,34};
CatmullRom(80) = {18,23,28,33};
CatmullRom(81) = {36,41,46,51};
CatmullRom(82) = {35,40,45,50};
CatmullRom(83) = {34,39,44,49};
CatmullRom(84) = {33,38,43,48};
CatmullRom(85) = {51,56,61,66};
CatmullRom(86) = {50,55,60,65};
CatmullRom(87) = {49,54,59,64};
CatmullRom(88) = {48,53,58,63};
CatmullRom(89) = {66,71,76,81};
CatmullRom(90) = {65,70,75,80};
CatmullRom(91) = {64,69,74,79};
CatmullRom(92) = {63,68,73,78};
CatmullRom(93) = {81,86,91,96};
CatmullRom(94) = {80,85,90,95};
CatmullRom(95) = {79,84,89,94};
CatmullRom(96) = {78,83,88,93};
CatmullRom(97) = {93,98,103,108};
CatmullRom(98) = {94,99,104,109};
CatmullRom(99) = {95,100,105,110};
CatmullRom(100) = {96,101,106,111};
CatmullRom(101) = {108,113,117};
CatmullRom(102) = {110,115,117};
CatmullRom(103) = {111,116,117};

// surfaces
Line Loop(104) = {70,-6,-71};
Ruled Surface(105) = {104};
Line Loop(106) = {71,-5,-4,-72};
Ruled Surface(107) = {106};
Line Loop(108) = {6,73,-12,-74};
Ruled Surface(109) = {108};
Line Loop(110) = {5,74,-11,-75};
Ruled Surface(111) = {110};
Line Loop(112) = {75,-10,-76,4};
Ruled Surface(113) = {112};
Line Loop(114) = {12,77,-21,-78};
Ruled Surface(115) = {114};
Line Loop(116) = {11,78,-20,-79};
Ruled Surface(117) = {116};
Line Loop(118) = {10,79,-19,-80};
Ruled Surface(119) = {118};
Line Loop(120) = {21,81,-30,-82};
Ruled Surface(121) = {120};
Line Loop(122) = {20,82,-29,-83};
Ruled Surface(123) = {122};
Line Loop(124) = {19,83,-28,-84};
Ruled Surface(125) = {124};
Line Loop(126) = {30,85,-39,-86};
Ruled Surface(127) = {126};
Line Loop(128) = {29,86,-38,-87};
Ruled Surface(129) = {128};
Line Loop(130) = {28,87,-37,-88};
Ruled Surface(131) = {130};
Line Loop(132) = {39,89,-48,-90};
Ruled Surface(133) = {132};
Line Loop(134) = {38,90,-47,-91};
Ruled Surface(135) = {134};
Line Loop(136) = {37,91,-46,-92};
Ruled Surface(137) = {136};
Line Loop(138) = {48,93,-57,-94};
Ruled Surface(139) = {138};
Line Loop(140) = {47,94,-56,-95};
Ruled Surface(141) = {140};
Line Loop(142) = {46,95,-55,-96};
Ruled Surface(143) = {142};
Line Loop(144) = {57,100,-66,-99};
Ruled Surface(145) = {144};
Line Loop(146) = {56,99,-65,-98};
Ruled Surface(147) = {146};
Line Loop(148) = {55,98,-64,-97};
Ruled Surface(149) = {148};
Line Loop(150) = {101,-102,-65,-64};
Ruled Surface(151) = {150};
Line Loop(152) = {103,-102,66};
Ruled Surface(153) = {152};
Symmetry {0,1,0,0} {
  Duplicata { Surface{105,107,111,109,113,119,117,115,121,123,125,127,129,131,133,135,137,139,141,143,145,147,149,153,151}; }
}
Symmetry {0,0,1,0} {
  Duplicata { Surface{154,158,105,107,111,113,109,168,173,163,119,117,115,188,183,178,203,198,193,121,123,125,213,208,127,129,131,218,233,228,223,133,135,137,143,141,139,238,243,248,263,258,253,145,147,149,272,268,153,151}; }
}
