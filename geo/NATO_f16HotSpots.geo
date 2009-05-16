lc = 0.01;
R1 = 0.201;
r1 = 0.0135;

// the base cylinder extrusion heights
h1 = 0.032;
h2 = h1 + 0.0366/2.0;
h3 = h1 + 0.0366;
h4 = 0.102;

Point(newp) = {0,0,0,lc};
Point(newp) = {0,-R1,0,lc};
Point(newp) = {0,-R1+0.058,0,lc};
Point(newp) = {-0.14,-R1+0.058,0,lc};
Point(newp) = {0,R1,0,lc};
Point(newp) = {-R1,0,0,lc};
Point(newp) = {0.0,-r1,0,lc};
Point(newp) = {r1,0,0,lc};
Point(newp) = {0,r1,0,lc};
Point(newp) = {-r1,0,0,lc};

// the tip
Point(newp) = {0.151,0,0,lc};
TipAngle = 51.5 / 180.0 * Pi;
LengthSideTip = 0.074;
Point(newp) = {0.151 - LengthSideTip * Cos(TipAngle/2.0), LengthSideTip * Sin(TipAngle/2.0),0,lc};
Point(newp) = {0.151 - LengthSideTip * Cos(TipAngle/2.0), -LengthSideTip * Sin(TipAngle/2.0),0,lc};
Point(newp) = {0.14,-R1+0.058,0,lc};
Point(newp) = {0.14,R1-0.058,0,lc};

// same points as above, but extruded following h2
Point(newp) = {-0.14,-R1+0.058,h2,lc};
Point(newp) = {-R1,0,h2,lc};
Point(newp) = {0,R1,h2,lc};
Point(newp) = {0.14,R1-0.058,h2,lc};

// same points as above, but extruded following h4
Point(newp) = {0,-R1,h4,lc};
Point(newp) = {0,-R1+0.058,h4,lc};
Point(newp) = {-0.14,-R1+0.058,h4,lc};
Point(newp) = {0,R1,h4,lc};
Point(newp) = {-R1,0,h4,lc};

Point(newp) = {0.151,0,h4,lc};
Point(newp) = {0.151 - LengthSideTip * Cos(TipAngle/2.0), LengthSideTip * Sin(TipAngle/2.0),h4,lc};
Point(newp) = {0.151 - LengthSideTip * Cos(TipAngle/2.0), -LengthSideTip * Sin(TipAngle/2.0),h4,lc};
Point(newp) = {0.14,-R1+0.058,h4,lc};
Point(newp) = {0.14,R1-0.058,h4,lc};


// the points for the 2 intakes on the lower cylinder
AngleBetweenIntakesEdges = 0.174/R1; // [rad]
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),0.0,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),0.0,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h4,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h4,lc};

AngleBetweenIntakesEdges = (0.174+0.0366)/R1; // [rad]
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),0.0,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),0.0,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h1,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h1,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h3,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h3,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h4,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h4,lc};

// the intakes
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0), (R1-0.0834) * Sin(AngleBetweenIntakesEdges/2.0),h1,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0), (R1-0.0834) * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0) - 0.0366/2.0 * Sin(AngleBetweenIntakesEdges/2.0), (R1-0.0834) * Sin(AngleBetweenIntakesEdges/2.0) - 0.0366/2.0 * Cos(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0) + 0.0366/2.0 * Sin(AngleBetweenIntakesEdges/2.0), (R1-0.0834) * Sin(AngleBetweenIntakesEdges/2.0) + 0.0366/2.0 * Cos(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0), (R1-0.0834) * Sin(AngleBetweenIntakesEdges/2.0),h3,lc};


Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0), (-R1+0.0834) * Sin(AngleBetweenIntakesEdges/2.0),h1,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0), (-R1+0.0834) * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0) + 0.0366/2.0 * Sin(AngleBetweenIntakesEdges/2.0), (-R1+0.0834) * Sin(AngleBetweenIntakesEdges/2.0) - 0.0366/2.0 * Cos(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0) - 0.0366/2.0 * Sin(AngleBetweenIntakesEdges/2.0), (-R1+0.0834) * Sin(AngleBetweenIntakesEdges/2.0) + 0.0366/2.0 * Cos(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {(-R1+0.0834) * Cos(AngleBetweenIntakesEdges/2.0), (-R1+0.0834) * Sin(AngleBetweenIntakesEdges/2.0),h3,lc};


AngleBetweenIntakesEdges = (0.174 + 2.0 * 0.0366)/R1; // [rad]
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),0.0,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),0.0,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h2,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), -R1 * Sin(AngleBetweenIntakesEdges/2.0),h4,lc};
Point(newp) = {-R1 * Cos(AngleBetweenIntakesEdges/2.0), R1 * Sin(AngleBetweenIntakesEdges/2.0),h4,lc};

// add the centres as we go up...
Point(newp) = {0.0, 0.0, h2,lc};
Point(newp) = {0.0, 0.0, h4,lc};

// the upper cylinder
R2 = 0.1; // 9.95 cm really :)
Point(newp) = {0.0506,0,h4,lc};
Point(newp) = {0.0506-R2,0,h4,lc};
Point(newp) = {0.0506-R2,-R2,h4,lc};
Point(newp) = {0.0506-R2,R2,h4,lc};
Point(newp) = {0.0506-0.1425,0,h4,lc};
Point(newp) = {0.0506-0.1425,-0.05,h4,lc};
Point(newp) = {0.0506-0.1425,0.05,h4,lc};

cos = 0.05/R2;
sin = Sqrt(1.0 - cos^2);
Point(newp) = {0.0506 -R2 - sin * R2,0.05,h4,lc};
Point(newp) = {0.0506 -R2 - sin * R2,-0.05,h4,lc};

// the curves for the hollow triangle in the upper cylinder
cos = 0.03/R2;
sin = Sqrt(1.0 - cos^2);
Point(newp) = {0.0506 -R2 - cos * R2, sin * R2,h4,lc};
Point(newp) = {0.0506 -R2 + cos * R2, sin * R2,h4,lc};



Line(1) = {2,3};
Line(3) = {14,13};
Line(4) = {13,11};
Line(5) = {11,12};
Line(6) = {12,15};
Line(7) = {15,19};
Line(8) = {19,29};
Line(9) = {29,26};
Line(10) = {26,25};
Line(11) = {25,27};
Line(12) = {27,28};
Line(14) = {21,20};
Line(15) = {20,2};
Line(16) = {3,21};
Line(17) = {28,14};
Line(18) = {27,13};
Line(19) = {25,11};
Line(20) = {12,26};


Line(21) = {3,4};
Line(22) = {22,21};
Line(23) = {22,16};
Line(24) = {16,4};
Line(25) = {60,58};
Line(26) = {58,56};
Line(27) = {36,38};
Line(28) = {42,44};
Line(29) = {34,32};
Line(30) = {32,30};
Line(31) = {24,17};
Line(32) = {17,6};
Line(33) = {35,33};
Line(34) = {33,31};
Line(35) = {45,43};
Line(36) = {39,37};
Line(37) = {61,59};
Line(38) = {59,57};
Line(39) = {23,18};
Line(40) = {18,5};
Line(41) = {53,58};
Line(42) = {38,51};
Line(43) = {54,32};
Line(44) = {55,42};
Line(45) = {50,43};
Line(46) = {48,33};
Line(47) = {46,39};
Line(48) = {49,59};
Line(49) = {8,13};
Line(50) = {8,12};
Line(51) = {9,5};
Line(52) = {7,3};
Line(53) = {10,6};

Line(54) = {69,70};
Line(55) = {69,72};
Line(56) = {70,71};
Circle(57) = {7,1,8};
Circle(58) = {8,1,9};

Circle(59) = {9,1,10};
Circle(60) = {10,1,7};
Circle(61) = {2,1,14};
Circle(62) = {15,1,5};
Circle(63) = {5,1,57};
Circle(64) = {57,1,37};
Circle(65) = {37,1,31};
Circle(66) = {31,1,6};
Circle(67) = {6,1,30};
Circle(68) = {30,1,36};
Circle(69) = {36,1,56};
Circle(70) = {56,1,4};
Circle(71) = {20,63,28};
Circle(72) = {29,63,23};
Circle(73) = {23,63,61};
Circle(74) = {61,63,45};
Circle(75) = {45,63,35};
Circle(76) = {35,63,24};
Circle(77) = {24,63,34};
Circle(78) = {34,63,44};
Circle(79) = {44,63,60};
Circle(80) = {60,63,22};
Circle(81) = {42,40,32};
Circle(82) = {32,40,38};
Circle(83) = {38,40,58};
Circle(84) = {58,40,42};
Circle(85) = {43,41,59};
Circle(86) = {59,41,39};
Circle(87) = {39,41,33};
Circle(88) = {33,41,43};
Circle(89) = {50,47,49};
Circle(90) = {49,47,46};
Circle(91) = {46,47,48};
Circle(92) = {48,47,50};
Circle(93) = {55,52,54};
Circle(94) = {54,52,51};
Circle(95) = {51,52,53};
Circle(96) = {53,52,55};

Circle(97) = {19,62,18};
Circle(98) = {18,62,59};
Circle(99) = {33,62,17};
Circle(100) = {17,62,32};
Circle(101) = {58,62,16};

// circles for the upper cylinder


Circle(102) = {64,65,66};
Circle(103) = {66,65,72};
Circle(104) = {71,65,73};
Circle(105) = {73,65,67};
Circle(106) = {67,65,74};
Circle(107) = {74,65,64};
Line Loop(108) = {22,-16,21,-24,-23};
Plane Surface(109) = {108};
Line Loop(110) = {16,14,15,1};
Plane Surface(111) = {110};
Line Loop(112) = {12,17,3,-18};
Plane Surface(113) = {112};
Line Loop(114) = {11,18,4,-19};
Plane Surface(115) = {114};
Line Loop(116) = {10,19,5,20};
Plane Surface(117) = {116};
Line Loop(118) = {6,7,8,9,-20};
Plane Surface(119) = {118};
Line Loop(120) = {49,4,5,-50};
Plane Surface(121) = {120};
Line Loop(122) = {52,-1,61,3,-49,-57};
Plane Surface(123) = {122};
Line Loop(124) = {51,-62,-6,-50,58};
Plane Surface(125) = {124};
Line Loop(126) = {53,-66,-65,-64,-63,-51,59};
Plane Surface(127) = {126};
Line Loop(128) = {53,67,68,69,70,-21,-52,-60};
Plane Surface(129) = {128};
Line(130) = {66,21};
Line(131) = {67,23};

Line(132) = {72,34};
Line(133) = {71,35};
Line(134) = {64,27};
Line(135) = {64,26};
Line Loop(136) = {130,14,71,-12,-134,102};
Plane Surface(137) = {136};
Line Loop(138) = {131,-72,9,-135,-107,-106};
Plane Surface(139) = {138};
Line Loop(140) = {131,73,74,75,-133,104,105};
Plane Surface(141) = {140};
Line Loop(142) = {22,-130,103,132,78,79,80};
Plane Surface(143) = {142};
Line Loop(144) = {132,-77,-76,-133,-56,-54,55};
Plane Surface(145) = {144};
Line Loop(146) = {89,90,91,92};
Plane Surface(147) = {146};
Line Loop(148) = {93,94,95,96};
Plane Surface(149) = {148};
Line Loop(150) = {72,39,-97,8};
Ruled Surface(151) = {150};
Line Loop(152) = {97,40,-62,7};
Ruled Surface(153) = {152};
Line Loop(154) = {73,37,-98,-39};
Ruled Surface(155) = {154};
Line Loop(156) = {98,38,-63,-40};
Ruled Surface(157) = {156};
Line Loop(158) = {76,31,-99,-33};
Ruled Surface(159) = {158};
Line Loop(160) = {99,32,-66,-34};
Ruled Surface(161) = {160};
Line Loop(162) = {77,29,-100,-31};
Ruled Surface(163) = {162};
Line Loop(164) = {100,30,-67,-32};
Ruled Surface(165) = {164};
Line Loop(166) = {80,23,-101,-25};
Ruled Surface(167) = {166};
Line Loop(168) = {101,24,-70,-26};
Ruled Surface(169) = {168};
Line Loop(170) = {78,-28,81,-29};
Ruled Surface(171) = {170};
Line Loop(172) = {79,25,84,28};
Ruled Surface(173) = {172};
Line Loop(174) = {82,-27,-68,-30};
Ruled Surface(175) = {174};
Line Loop(176) = {27,83,26,-69};
Ruled Surface(177) = {176};
Line Loop(178) = {65,-34,-87,36};
Ruled Surface(179) = {178};
Line Loop(180) = {36,-64,-38,86};
Ruled Surface(181) = {180};
Line Loop(182) = {74,35,85,-37};
Ruled Surface(183) = {182};
Line Loop(184) = {75,33,88,-35};
Ruled Surface(185) = {184};
Line Loop(186) = {45,85,-48,-89};
Ruled Surface(187) = {186};
Line Loop(188) = {45,-88,-46,92};
Ruled Surface(189) = {188};
Line Loop(190) = {46,-87,-47,91};
Ruled Surface(191) = {190};
Line Loop(192) = {48,86,-47,-90};
Ruled Surface(193) = {192};
Line Loop(194) = {44,81,-43,-93};
Ruled Surface(195) = {194};
Line Loop(196) = {96,44,-84,-41};
Ruled Surface(197) = {196};
Line Loop(198) = {43,82,42,-94};
Ruled Surface(199) = {198};
Line Loop(200) = {95,41,-83,42};
Ruled Surface(201) = {200};
Line Loop(202) = {135,10,11,-134};
Plane Surface(203) = {202};
Line Loop(204) = {71,17,-61,-15};
Ruled Surface(205) = {204};
